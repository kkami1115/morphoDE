#!/usr/bin/env python3
"""
morphoDE dataset downloader — reproducible fetch of the 11 analysis datasets.

Reads ../data/catalog.csv and fetches one dataset (or all) from its primary
public source. Downloads resume on interruption and are verified: the first
successful fetch of a file records its sha256 in a per-dataset lockfile
(data/<id>/CHECKSUMS.sha256); re-running verifies existing files against that
lockfile and skips ones that already match. Where the upstream repository
publishes its own checksum (Zenodo, Figshare), the fetched file is additionally
checked against it.

This does NOT redistribute data. It fetches from the original repositories, so
the exact bytes are whatever those repositories currently serve. The lockfile
makes YOUR download reproducible and lets you detect a corrupted or changed file
on re-run; it is written on first run, not shipped pre-populated.

Sources handled automatically (catalog column `source`):
  geo    - NCBI GEO series supplementary files      (source_arg = GSExxxxxx)
  zenodo - Zenodo record                            (source_arg = record id)
  figshare - Figshare article                       (source_arg = article id)
  cngb   - CNGB / STOmicsDB collection              (source_arg = STDS/STT id,
                                                      source_arg2 = portal slug)

Datasets with `source = manual` (or `access = manual`) cannot be fetched by a
single command — controlled-access, or a raw-data project with no stable direct
link. For those the script prints the accession, DOI and portal URL and stops,
so you obtain them through the proper channel. This is honest by design: it will
not silently produce a partial mirror and call it complete.

Usage:
  python3 download.py --list                 # show catalog and access class
  python3 download.py <dataset_id>           # fetch one dataset
  python3 download.py --all                  # fetch every 'auto' dataset
  python3 download.py <dataset_id> --verify  # re-check local files vs lockfile
  python3 download.py <dataset_id> --dry-run # resolve & list URLs, fetch nothing

Options:
  --data-dir PATH   destination root (default: ../data, or $MORPHODE_DATA_DIR)
  --workers N       parallel file downloads (default 4)
  --compress        OPTIONAL: zstd-compress *.h5ad after download to save disk
                    (needs `zstd`). Off by default; the datasets download as
                    plain .h5ad and the analysis reads plain .h5ad. This flag is
                    only a local disk-space convenience and is not required for
                    reproduction.

Requires: python3, curl. Optional: zstd (for --compress).
"""
import argparse, csv, hashlib, os, re, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CATALOG = os.path.join(REPO, "data", "catalog.csv")
EXT = r'\.(?:h5ad|h5|gef|gem|loom|csv|tsv|txt|gz|tar|tar\.gz|zip|rds|mtx|json)'
_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print("[%s] %s" % (datetime.now().strftime("%H:%M:%S"), msg), flush=True)


def load_catalog():
    with open(CATALOG, newline="") as fh:
        return {r["dataset_id"]: r for r in csv.DictReader(fh)}


# ---- source resolvers: (dataset row) -> list of (relpath, url, upstream_md5) ----

def curl_text(url, tries=4):
    for _ in range(tries):
        p = subprocess.run(["curl", "-sL", "--max-time", "90", "--retry", "3", url],
                            capture_output=True, text=True)
        if p.returncode == 0 and p.stdout:
            return p.stdout
    return None


def curl_json(url):
    import json
    t = curl_text(url)
    if not t:
        return None
    try:
        return json.loads(t)
    except Exception:
        return None


def resolve_geo(row):
    gse = row["source_arg"]
    base = "https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/" % (gse[:-3], gse)
    html = curl_text(base)
    if html is None:
        sys.exit("GEO listing failed for %s" % gse)
    names = sorted(set(re.findall(r'href="([^"?/][^"]*%s)"' % EXT, html, re.I)))
    return [(n, base + n, "") for n in names]


def resolve_zenodo(row):
    rid = row["source_arg"]
    meta = curl_json("https://zenodo.org/api/records/%s" % rid)
    if not meta or "files" not in meta:
        sys.exit("Zenodo metadata failed for record %s" % rid)
    out = []
    for f in meta["files"]:
        k = f["key"]
        out.append((k, "https://zenodo.org/api/records/%s/files/%s/content" % (rid, k),
                    (f.get("checksum") or "").replace("md5:", "")))
    return out


def resolve_figshare(row):
    aid = row["source_arg"]
    meta = curl_json("https://api.figshare.com/v2/articles/%s" % aid)
    if not meta or "files" not in meta:
        sys.exit("Figshare metadata failed for article %s" % aid)
    return [(f["name"], f["download_url"],
             f.get("computed_md5") or f.get("supplied_md5") or "") for f in meta["files"]]


def _walk_cngb_index(base, depth=0, max_depth=4):
    """Recursively enumerate file URLs under a CNGB SciRAID index directory.
    The FTP-over-HTTPS index is a plain autoindex: files (matching EXT) are
    collected; subdirectories (href ending in '/') are descended, since some
    collections nest their .h5ad under a subdir (e.g. STDS0000104/stomics/)."""
    idx = curl_text(base)
    if not idx or depth > max_depth:
        return []
    files, subdirs = [], []
    for href in re.findall(r'href="([^"?]+)"', idx, re.I):
        if href.startswith("/") or href.startswith("http") or href.startswith(".."):
            continue  # site chrome / parent link, not a tree entry
        full = base + href
        if href.endswith("/"):
            subdirs.append(full)
        elif re.search(EXT + r"$", href, re.I):
            files.append(full)
    for sub in subdirs:
        files.extend(_walk_cngb_index(sub, depth + 1, max_depth))
    return sorted(set(files))


def resolve_cngb(row):
    """CNGB / STOmicsDB. Prefer the server-rendered portal download page
    (db.cngb.org/stomics/<slug>/download/) which stays up when the FTP index
    refuses connections; fall back to enumerating the SciRAID FTP tree by STDS id.
    Files are served from ftp.cngb.org; we rewrite to ftp:// which is more reliable."""
    stds, slug = row["source_arg"], row.get("source_arg2", "")
    file_re = (r'https://ftp\.cngb\.org/pub/SciRAID/stomics/'
               r'(?:STDS[0-9]{7}|STT[0-9]{7})/[^"\'\s<>]+%s' % EXT)
    urls = []
    if slug:
        page = curl_text("https://db.cngb.org/stomics/%s/download/" % slug)
        if page:
            urls = sorted(set(re.findall(file_re, page)))
    if not urls:  # fall back to a recursive walk of the SciRAID index tree
        base = "https://ftp.cngb.org/pub/SciRAID/stomics/%s/" % stds
        urls = _walk_cngb_index(base)
    if not urls:
        sys.exit("CNGB enumeration failed for %s (slug=%s). Fetch manually from %s"
                 % (stds, slug, row["url"]))
    root = "/%s/" % stds
    out = []
    for u in urls:
        rel = u.split(root, 1)[1] if root in u else u.rsplit("/", 1)[-1]
        # keep the canonical https URL; curl_download tries ftp:// as a fallback
        # (ftp.cngb.org's https and ftp fronts each flap independently by network)
        out.append((rel, u if u.startswith("http") else "https:" + u.split(":", 1)[1], ""))
    return out


RESOLVERS = {"geo": resolve_geo, "zenodo": resolve_zenodo,
             "figshare": resolve_figshare, "cngb": resolve_cngb}


# ---- checksum / lockfile ----

def sha256_file(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_file(path, buf=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def load_lock(dst_dir):
    lock = {}
    p = os.path.join(dst_dir, "CHECKSUMS.sha256")
    if os.path.exists(p):
        for line in open(p):
            line = line.strip()
            if line and "  " in line:
                digest, rel = line.split("  ", 1)
                lock[rel] = digest
    return lock


def write_lock(dst_dir, lock):
    p = os.path.join(dst_dir, "CHECKSUMS.sha256")
    with open(p, "w") as fh:
        for rel in sorted(lock):
            fh.write("%s  %s\n" % (lock[rel], rel))


# ---- fetch ----

def _curl_once(url, dest):
    return subprocess.run(
        ["curl", "-sL", "--fail", "-C", "-", "--retry", "10", "--retry-delay", "5",
         "--retry-all-errors", "--connect-timeout", "30",
         "--speed-limit", "1000", "--speed-time", "120", "-o", dest, url]
    ).returncode == 0


def curl_download(url, dest):
    """Download `url` to `dest`. For CNGB (ftp.cngb.org) the https and native-ftp
    fronts each flap independently depending on the network, so try both schemes:
    whichever the caller gave first, then the other."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    candidates = [url]
    if "ftp.cngb.org" in url:
        if url.startswith("https://"):
            candidates.append("ftp://" + url[len("https://"):])
        elif url.startswith("ftp://"):
            candidates.append("https://" + url[len("ftp://"):])
    for u in candidates:
        if _curl_once(u, dest) and not _is_html_stub(dest):
            return True
    return False


def _is_html_stub(path):
    """A flapping CNGB front sometimes returns a tiny HTML redirect page with a
    2xx code instead of the file. Reject it so the other scheme is tried."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(512).lstrip().lower()
        return head.startswith(b"<!doctype html") or head.startswith(b"<html")
    except OSError:
        return True


def fetch_one(task, dst_dir, lock, verify_only):
    rel, url, up_md5 = task
    dest = os.path.join(dst_dir, rel)
    have = os.path.exists(dest)
    # already verified against lockfile?
    if have and rel in lock:
        if sha256_file(dest) == lock[rel]:
            return (rel, "OK_CACHED")
        return (rel, "CHANGED" if verify_only else "REDOWNLOAD")
    if verify_only:
        return (rel, "MISSING" if not have else "UNTRACKED")
    if not curl_download(url, dest):
        return (rel, "DL_FAIL")
    if up_md5 and md5_file(dest) != up_md5:
        return (rel, "MD5_MISMATCH")
    lock[rel] = sha256_file(dest)
    return (rel, "OK")


def maybe_compress(dst_dir):
    if not any(subprocess.run(["which", "zstd"], capture_output=True).stdout):
        log("zstd not found; skipping --compress"); return
    for root, _, files in os.walk(dst_dir):
        for f in files:
            if f.endswith(".h5ad"):
                p = os.path.join(root, f)
                subprocess.run(["zstd", "-19", "--long=27", "-T0", "-q", "--check", "--rm", p])


def do_dataset(row, data_dir, workers, verify_only, dry_run, compress):
    did = row["dataset_id"]
    if row["access"].startswith("manual") or row["source"] == "manual":
        log("MANUAL: %s is not one-command downloadable." % did)
        print("  reason : %s" % row["notes"])
        print("  accession: %s" % row["accession"])
        print("  DOI    : https://doi.org/%s" % row["doi"])
        print("  portal : %s" % row["url"])
        print("  license: %s" % row["license"])
        return
    resolver = RESOLVERS.get(row["source"])
    if resolver is None:
        sys.exit("no resolver for source=%s (%s)" % (row["source"], did))
    log("resolving %s (%s %s)" % (did, row["source"], row["source_arg"]))
    tasks = resolver(row)
    log("  %d file(s) at source" % len(tasks))
    if dry_run:
        for rel, url, md5 in tasks:
            print("  %s\n      %s%s" % (rel, url, "  md5=" + md5 if md5 else ""))
        return
    dst_dir = os.path.join(data_dir, did)
    os.makedirs(dst_dir, exist_ok=True)
    lock = load_lock(dst_dir)
    results = []
    lock_lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_one, t, dst_dir, lock, verify_only): t[0] for t in tasks}
        for fu in as_completed(futs):
            rel, status = fu.result()
            results.append((rel, status))
            log("  [%s] %s" % (status, rel))
            # Persist the lockfile incrementally as each file lands, so an
            # interrupted run keeps its verified progress (datasets can be 100s
            # of files) instead of re-downloading everything on the next run.
            if not verify_only and status == "OK":
                with lock_lock:
                    write_lock(dst_dir, lock)
    if not verify_only:
        write_lock(dst_dir, lock)
        if compress:
            maybe_compress(dst_dir)
    from collections import Counter
    log("%s summary: %s" % (did, dict(Counter(s for _, s in results))))


def main():
    ap = argparse.ArgumentParser(description="morphoDE reproducible dataset downloader")
    ap.add_argument("dataset_id", nargs="?", help="catalog dataset_id, or use --all")
    ap.add_argument("--all", action="store_true", help="fetch every 'auto' dataset")
    ap.add_argument("--list", action="store_true", help="print catalog and exit")
    ap.add_argument("--verify", action="store_true", help="check local files vs lockfile; download nothing")
    ap.add_argument("--dry-run", action="store_true", help="resolve URLs and list; download nothing")
    ap.add_argument("--compress", action="store_true",
                    help="OPTIONAL local disk-space convenience: zstd-compress *.h5ad after "
                         "download. Off by default; not required for reproduction.")
    ap.add_argument("--data-dir", default=os.environ.get("MORPHODE_DATA_DIR", os.path.join(REPO, "data")))
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    cat = load_catalog()
    if args.list or (not args.dataset_id and not args.all):
        print("%-34s %-12s %-8s %s" % ("dataset_id", "source", "access", "accession"))
        print("-" * 90)
        for did, r in cat.items():
            print("%-34s %-12s %-8s %s" % (did, r["source"], r["access"], r["accession"]))
        if not args.list:
            print("\nUsage: python3 download.py <dataset_id> | --all | --list")
        return

    targets = list(cat) if args.all else [args.dataset_id]
    for did in targets:
        if did not in cat:
            sys.exit("unknown dataset_id: %s (see --list)" % did)
        do_dataset(cat[did], args.data_dir, args.workers, args.verify, args.dry_run, args.compress)


if __name__ == "__main__":
    main()

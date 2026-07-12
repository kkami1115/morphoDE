# morphoDE — 3-minute demo video script

**Format:** hybrid — figure/narration for the finding, live Claude Science screen capture for the method and reproducibility.
**Target length:** 180 s (hard cap 3:00).
**Narration language:** English (judges are English-speaking). A Japanese reference translation is given under each block for recording; use whichever you record in — keep on-screen text English either way.
**Aspect / capture:** 1920×1080. Record the Claude Science session at the same resolution so text stays legible when cut in.

Every headline number below is verified against `results/bd_final_headline.csv` (126 interface–strata): signal in **125/126**, **97** at ≥6/12, **16** perfect 12/12, **1** zero, signal in **11/11** datasets.

---

## Shot list

### 1 · Hook — the question (0:00–0:18, ~18 s)
**On screen:** Black title card → fade to the `crossspecies_tree_of_life.png` headline figure, still, not yet explained. Title text: *"Does the shape of a tissue interface change what its cells express?"*

**Narration (EN):**
> "Tissues are full of interfaces — where a tumor meets healthy tissue, where one cell type meets another. Some parts of that boundary bulge outward; others cave inward. A simple question: do the cells on a bulge express different genes than the cells in a dent? We tested it across the tree of life."

**日本語(参考):** 組織は界面だらけです——腫瘍と正常組織の境目、細胞種どうしの境目。その境界には外へ出っ張る部分と内へ引っ込む部分がある。素朴な問い:出っ張りの細胞は、引っ込みの細胞と違う遺伝子を出しているのか? これを生命の樹全体で検証しました。

---

### 2 · The finding (0:18–0:52, ~34 s)
**On screen:** Animate/zoom the headline figure. Highlight (a) the phylogenetic row order spanning human→Arabidopsis, (b) dots pushing to the right, (c) the many interfaces reaching 12/12 across unrelated species.

**Narration (EN):**
> "This is the answer. Eleven 3-D spatial-transcriptomics datasets, a hundred and twenty-six interface–strata, from a human lymph-node tumor to an Arabidopsis leaf. Each dot is one interface in one specimen or one time point; the further right, the more genes whose bulge-versus-dent difference survives a strict null. A hundred and twenty-five of a hundred and twenty-six carry signal — every dataset, every branch of the tree. Sixteen reach a perfect twelve out of twelve, in species that share no common ancestry. And because we tested every developmental stage and every regeneration time point separately, this isn't one lucky snapshot — it holds across conditions. The effect is not one tissue's quirk. It recurs across the tree of life."

**日本語(参考):** これが答えです。11の3次元空間トランスクリプトーム・データセット、126の界面-層、ヒトのリンパ節腫瘍からシロイヌナズナの葉まで。各点は1検体・1時点における1つの界面で、右にあるほど、出っ張り対引っ込みの発現差が厳しいヌルを生き残った遺伝子が多い。126中125が信号あり——全データセット、生命樹の全枝で。16が12分の12満点で、系統的に無関係な種にまたがる。しかも各発生段階・各再生時点を別々に検証したので、これは偶然の1スナップショットではなく条件を越えて成立する。特定組織の癖ではなく、生命の樹を越えて再現する現象です。

---

### 3 · How it's measured — and why it's honest (0:52–1:28, ~36 s)
**On screen:** Cut to `tumor_interface_proof.png`. Walk panels in order: A (open interface sheet, two angles) → B (93% within cell type) → C (within-tumor DE, bulge/convex=immunoglobulin/cholesterol, dent/concave=hypoxia/invasion) → D (geometry-breaking null) → E (power curve). Keep each on screen ~7 s.

**Narration (EN):**
> "Here's how, on the human tumor. The interface is one open curved sheet — the real tumor-to-stroma boundary, not a closed blob — and we score each patch by how far it bulges out or caves in, using its normal, no curvature and no arbitrary threshold. Ninety-three percent of the difference is within a single tumor cell type, so this is a change in cell state, not a reshuffle of cell types. On the convex bulges: immunoglobulin, antigen presentation, cholesterol synthesis. In the concave dents: hypoxia, invasion, keratinization — PTHLH, KRT17. The honest part is the null: we shuffle the geometry two hundred times, keeping the biology fixed. Real signal survives; a coincidence of shape would not. The test is calibrated at five percent and powered to a log-fold of two-hundredths."

**日本語(参考):** ヒト腫瘍でのやり方です。界面は1枚の開いた曲面——閉じた塊ではなく、本物の腫瘍-間質境界——で、各パッチが法線方向にどれだけ出っ張るか引っ込むかで採点する。曲率も恣意的な閾値も使わない。差の93%が単一腫瘍細胞種の中にあるので、これは細胞種の入れ替えでなく細胞状態の変化。凸の出っ張りには免疫グロブリン・抗原提示・コレステロール合成、凹の引っ込みには低酸素・浸潤・角化(PTHLH, KRT17)。誠実さの肝はヌル:幾何だけを200回壊し、生物学は固定したまま。本物の信号は生き残り、形の偶然なら生き残らない。検定は5%で較正され、log-fold 0.02の検出力があります。

---

### 4 · Claude Use — reframe, self-correct, scale, shown live (1:28–2:32, ~64 s)
**On screen:** Switch to a screen capture of the Claude Science session. Three beats, each visible on screen:
1. Scroll to where the naive correlation/SAR analysis was failing, then the message reframing it into the geometry-breaking null. *(caption overlay: "1 — Claude reframed a failing test into a geometry-breaking null")*
2. Scroll to the moment the reviewer flags the closed-blob interface, then Claude's quantitative check (bulge cells face ~4% non-tumor — empty space, not the boundary) and the rebuild on the open sheet. *(caption: "2 — caught its own structural flaw and re-ran all 126 interface–strata")*
3. Open the method `open_interface.py` / skill and run it on one dataset with a single call, showing the open sheet, DE and null tables appear. *(caption: "3 — packaged as a reusable method — one call on a new dataset")*

**Narration (EN):**
> "Claude didn't just run this — it shaped it, and it corrected itself. The first analysis was failing to separate real signal from tissue shape; Claude proposed the fix — break the geometry, keep the biology, so the null answers exactly the confound we worried about. Then a harder moment: the interface was being drawn as a closed blob when the tumor was never fully acquired. Instead of trusting the figure, Claude tested it — the bulges were facing empty, unsampled space, not the real boundary — and rebuilt the whole method on the genuine open surface, re-running all a hundred and twenty-six interfaces. The finding didn't just survive; it sharpened. And it's packaged: one call, a new dataset, and the surface, the null, and the differential expression come back. That's a workflow a reviewer can rerun, not just a figure to trust."

**日本語(参考):** Claudeはこれを実行しただけでなく、形にし、自己修正した。最初の解析は本物の信号と組織の形を分離できずにいた——Claudeが修正を提案:幾何を壊し生物学を保つ。そしてより難しい局面:腫瘍は完全には取得されていないのに界面が閉じた塊として描かれていた。図を信じる代わりにClaudeは検証し——出っ張りが実境界でなく未サンプルの空間に面していた——手法全体を真の開曲面に作り直し、126界面を全て再実行した。発見は生き残っただけでなく、鋭くなった。しかも再利用可能:1コール、新データセット、そして曲面・ヌル・発現差が返る。信じてもらう図ではなく、査読者が再実行できるワークフローです。

---

### 5 · Impact + close (2:32–3:00, ~28 s)
**On screen:** Split screen — headline figure on the left, the GitHub repo page (github.com/kkami1115/morphoDE) on the right. End on the repo URL + the one-line claim.

**Narration (EN):**
> "So interface shape and cell state are linked — and the link is general, from cancer to plants. In the tumor, the concave dents carry the hypoxia-and-invasion program and the convex bulges carry immune and cholesterol signatures — a spatial map of tumor cell state read straight off geometry. Everything — the method, the eleven datasets, the null, the figures — is open on GitHub. That's morphoDE: bulge versus dent, tested across the tree of life."

**日本語(参考):** つまり界面の形と細胞の状態は結びついており、その結びつきはがんから植物まで一般的です。腫瘍では、凹の引っ込みが低酸素・浸潤プログラムを、凸の出っ張りが免疫・コレステロールのシグネチャを担う——腫瘍細胞状態の空間地図が幾何から直接読める。手法、11データセット、ヌル、図——すべてGitHubで公開。これがmorphoDE:出っ張り対引っ込みを、生命の樹を越えて検証した仕事です。

---

## Production notes
- **Timing budget:** 18 + 34 + 36 + 64 + 28 = 180 s exactly. The live segment (§4) is the compressible one — if a beat runs long, trim beat 1's scroll, not the skill re-run (that's the Claude-Use payload).
- **Legibility:** the six-panel proof is dense; when a panel is on screen, crop/zoom to just that panel rather than showing the whole figure small.
- **Live re-run choice:** pick a *fast* dataset for the on-camera run (zebrafish Yolk Syncytial Layer or axolotl MSN at 10DPI — both small and both perfect 12/12, so the result table is clean and impressive). Avoid planarian (millions of cells — too slow on camera).
- **Captions:** burn the three §4 captions in as on-screen text; judges may watch muted.
- **What to say on-screen vs. narrate:** all figure titles are already English, so let them carry the labels; narration adds the interpretation, not the axis names.
- **Sign is easy to say wrong:** on camera, always pair the word with the geometry — "convex bulge = immune/cholesterol", "concave dent = hypoxia/invasion". The invasion program is on the DENT (concave) side; do not call the bulge the invasion signature.
- **Cold-open alternative:** if the hook feels slow, open on the 3-D open interface sheet rotating (`concept_3d_rotate.mp4`) for 3 s before the title card — motion holds attention.

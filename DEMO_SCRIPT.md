# morphoDE — 3-minute demo video script

**Format:** hybrid — figure/narration for the finding, live Claude Science screen capture for the method and reproducibility.
**Target length:** 180 s (hard cap 3:00).
**Narration language:** English (judges are English-speaking). A Japanese reference translation is given under each block for recording; use whichever you record in — keep on-screen text English either way.
**Aspect / capture:** 1920×1080. Record the Claude Science session at the same resolution so text stays legible when cut in.

Every headline number below is verified against `results/bd_all_meta.csv` (38 interfaces): signal in **35/38** interfaces, **19** at ≥6/12, **4** perfect 12/12, **0** in 3, signal in **11/11** datasets.

---

## Shot list

### 1 · Hook — the question (0:00–0:18, ~18 s)
**On screen:** Black title card → fade to the `crossspecies_tree_of_life.png` headline figure, still, not yet explained. Title text: *"Does the shape of a tissue interface change what its cells express?"*

**Narration (EN):**
> "Tissues are full of interfaces — where a tumor meets healthy tissue, where one cell type meets another. Some parts of that boundary bulge outward; others cave inward. A simple question: do the cells on a bulge express different genes than the cells in a dent? We tested it across the tree of life."

**日本語(参考):** 組織は界面だらけです——腫瘍と正常組織の境目、細胞種どうしの境目。その境界には外へ出っ張る部分と内へ引っ込む部分がある。素朴な問い:出っ張りの細胞は、引っ込みの細胞と違う遺伝子を出しているのか? これを生命の樹全体で検証しました。

---

### 2 · The finding (0:18–0:52, ~34 s)
**On screen:** Animate/zoom the headline figure. Highlight (a) the cladogram spanning human→Arabidopsis, (b) dots pushing to the right, (c) the four labelled 12/12 interfaces.

**Narration (EN):**
> "This is the answer. Eleven 3-D spatial-transcriptomics datasets, thirty-eight interfaces, from a human lymph-node tumor to the tip of an Arabidopsis root. Each dot is one interface; the further right, the more genes whose bulge-versus-dent difference survives a strict null. Signal shows up in every single dataset — thirty-five of thirty-eight interfaces. And four interfaces that share no common ancestry — a human tumor, an axolotl vascular layer, the zebrafish notochord, a planarian blastema — hit a perfect twelve out of twelve. The effect is not one tissue's quirk. It recurs across the tree of life."

**日本語(参考):** これが答えです。11の3次元空間トランスクリプトーム・データセット、38の界面、ヒトのリンパ節腫瘍からシロイヌナズナの根端まで。各点が1つの界面で、右にあるほど、出っ張り対引っ込みの発現差が厳しいヌルを生き残った遺伝子が多い。信号は全データセットで出た——38界面中35。しかも系統的に無関係な4界面(ヒト腫瘍・アホロートル血管層・ゼブラフィッシュ脊索・プラナリア再生芽)が12分の12満点。特定組織の癖ではなく、生命の樹を越えて再現する現象です。

---

### 3 · How it's measured — and why it's honest (0:52–1:28, ~36 s)
**On screen:** Cut to `tumor_interface_proof.png`. Walk panels in order: A (3-D mesh, two angles) → C (within-tumor DE, bulge=hypoxia/invasion, dent=immunoglobulin) → D (geometry-breaking null) → E (power curve). Keep each on screen ~8 s.

**Narration (EN):**
> "Here's how, on the human tumor. We reconstruct the interface as a 3-D surface and score each patch by how far it deviates from the tissue's convex hull — no curvature, no arbitrary threshold. Bulge cells run high in hypoxia and invasion genes — SLC2A1, EGLN3, PTHLH; dent cells run high in immunoglobulins. The honest part is the null: we shuffle the geometry — rotate the bulge/dent labels around each section — two hundred times, keeping the biology fixed. Real signal survives; a coincidence of shape would not. The test is calibrated at five percent and powered to detect a log-fold of two-hundredths."

**日本語(参考):** ヒト腫瘍でのやり方です。界面を3次元曲面として再構成し、各パッチが組織の凸包からどれだけ外れるかで採点する——曲率も恣意的な閾値も使わない。出っ張りの細胞は低酸素・浸潤の遺伝子(SLC2A1, EGLN3, PTHLH)が高く、引っ込みの細胞は免疫グロブリンが高い。誠実さの肝はヌル:幾何だけを壊す——各切片で出っ張り/引っ込みラベルを回転させる——のを200回、生物学は固定したまま。本物の信号は生き残り、形の偶然なら生き残らない。検定は5%で較正され、log-fold 0.02を検出する検出力があります。

---

### 4 · Claude Use — reframe + scale, shown live (1:28–2:32, ~64 s)
**On screen:** Switch to a screen capture of the Claude Science session. Three beats, each visible on screen:
1. Scroll to the moment the naive correlation/SAR analysis was failing, then the message where the analysis is reframed into the geometry-breaking null. *(caption overlay: "1 — Claude reframed a failing test into a geometry-breaking null")*
2. Show the parallel fan-out running the pipeline across all 11 datasets. *(caption: "2 — fanned the pipeline across 11 datasets in parallel")*
3. Open the published skill `interface-bulge-dent-de` and run it on one dataset with a single command, showing the meta/null/DE tables appear. *(caption: "3 — packaged as a reusable skill — one command on a new dataset")*

**Narration (EN):**
> "Claude didn't just run this — it shaped it. The first analysis, a naive correlation, was failing to separate real signal from tissue shape. Claude proposed the fix: break the geometry and keep the biology, so the null answers exactly the confound we were worried about. Then it scaled — running the same pipeline across all eleven datasets in parallel — and packaged the whole method as a reusable skill. Watch: one command, a new dataset, and the mesh, the null, and the differential expression come back. That's the workflow a reviewer can rerun, not just a figure to trust."

**日本語(参考):** Claudeはこれを実行しただけでなく、形にした。最初の解析(素朴な相関)は、本物の信号と組織の形を分離できずにいた。Claudeが修正を提案:幾何を壊し生物学を保つ——こうすればヌルが、まさに我々が懸念した交絡に答える。そして拡張——同じパイプラインを11データセットで並列実行——手法全体を再利用可能なskillにまとめた。ご覧の通り:1コマンド、新しいデータセット、そしてメッシュ・ヌル・発現差が返ってくる。信じてもらう図ではなく、査読者が再実行できるワークフローです。

---

### 5 · Impact + close (2:32–3:00, ~28 s)
**On screen:** Split screen — headline figure on the left, the GitHub repo page (github.com/kkami1115/morphoDE) on the right. End on the repo URL + the one-line claim.

**Narration (EN):**
> "So interface shape and cell state are linked — and the link is general, from cancer to plants. For cancer, the bulge program is exactly the invasion-and-hypoxia signature, read straight off tissue geometry. Everything — the skill, the eleven datasets, the null, all five figures — is open on GitHub. That's morphoDE: bulge versus dent, tested across the tree of life."

**日本語(参考):** つまり界面の形と細胞の状態は結びついており、その結びつきはがんから植物まで一般的です。がんでは、出っ張りのプログラムがまさに浸潤・低酸素のシグネチャで、組織の形から直接読める。skill、11データセット、ヌル、5枚の図——すべてGitHubで公開。これがmorphoDE:出っ張り対引っ込みを、生命の樹を越えて検証した仕事です。

---

## Production notes
- **Timing budget:** 18 + 34 + 36 + 64 + 28 = 180 s exactly. The live segment (§4) is the compressible one — if a beat runs long, trim beat 1's scroll, not the skill re-run (that's the Claude-Use payload).
- **Legibility:** the six-panel proof is dense; when a panel is on screen, crop/zoom to just that panel rather than showing the whole figure small.
- **Live re-run choice:** pick a *fast* dataset for the on-camera skill run (zebrafish notochord = 3,917 cells, or axolotl VLMC = 28,661 — both small and both perfect 12/12, so the result table is clean and impressive). Avoid planarian l4 (2.7 M cells — too slow on camera).
- **Captions:** burn the three §4 captions in as on-screen text; judges may watch muted.
- **What to say on-screen vs. narrate:** all figure titles are already English, so let them carry the labels; narration adds the interpretation, not the axis names.
- **Cold-open alternative:** if the hook feels slow, open on the 3-D tumor mesh rotating (panel A) for 3 s before the title card — motion holds attention.

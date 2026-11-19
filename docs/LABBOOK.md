# autogel – Lab book

Development record: the original idea, decisions with their reasons, and a dated log of changes.
New entries go at the **top** of the change log, and every change entry says *what*, *why* and
*how it was checked*.

---

## 1. Initial idea (2026-11-18)

Requirements as first stated:
- automatic band and lane detection and band intensities on agarose and SDS gels
- user input: the sizes of the marker
- the user can correct the automatic recognition (lanes, bands, …)
- output: table of intensities, log of what was done, intensity profile of each lane
- if possible build on ImageJ/Fiji functions, so Python
- open questions: does something comparable exist? what is missing?

Requirements added during the first session:
- ladder library with vendor ladders plus user-defined ones
- check the lane count with the user before detection runs
- documentation: a user manual and this lab book
- a git project
- light on RAM and GPU, so it runs on a budget PC or old laptop
- readable labels; colour-blind-safe colours
- process rule: lint, review and stress-test before every commit

## 2. Existing tools (survey, 2026-11-18)

| Tool | Notes |
|---|---|
| ImageJ/Fiji *Analyze › Gels* | manual lanes and peak areas; no automatic detection or sizing |
| GelJ (Fiji-based) | lane/band detection and sizing, aimed at DNA fingerprinting |
| GelAnalyzer (free, Windows) | automatic, closed source, can't be scripted |
| GelGenie (2024/25, open source) | deep-learning segmentation; limited sizing and editing workflow |
| Image Lab, ImageQuant TL, Image Studio, TotalLab | commercial, complete |

Gap: an open, scriptable, light-weight tool that is both automatic and correctable, with a
replayable log.

## 3. Design decisions

| # | Decision | Reason |
|---|---|---|
| D1 | Core in numpy/scipy/scikit-image, **no Fiji dependency** | ImageJ's gel tools are manual GUI steps; the useful algorithms (e.g. rolling ball) exist in Python. Avoids Java and the Fiji install. PyImageJ stays an option for ImageJ-identical results |
| D2 | GUI = **napari** plugin | Chosen by the user. Native Python, editable shape layers for lanes/bands |
| D3 | Orientation from **gradient energy** (bands, wells and tray edges all run across the lanes) | The first attempt, periodicity of the column profile, was fooled by the gel-tray edges in Test_1 and rotated the gel the wrong way |
| D4 | Lane finding on a **vertical top-hat** image (keeps short horizontal structures) | Removes tray edges, smears and illumination gradients, which created false lanes |
| D5 | Lane spacing from the **first major autocorrelation peak**, refined on a comb-grid score. Centres from the half-maximum plateau of each lane | A grid score alone favours 3×/5× spacings, which hit only the bright ladder lanes. A plain sum favours tiny spacings. A fixed grid drifts across 20 lanes |
| D6 | **Ask for or confirm the lane count** (GUI), `--lanes` (CLI) | The user's idea. With n known, the detector only picks which run of lanes is the gel, so tray edges and stickers are excluded |
| D7 | Band integration on a lane profile = **sum across 80 % of the lane width** per row, minus a 1-D rolling-ball baseline | Standard densitometry; fast. 1-D baseline instead of a 2-D rolling ball on 5 MP keeps it light |
| D8 | **Lane-coverage filter** (band must cover ≥ 40 % of the lane width) | Dust and scratches gave extra "bands" (e.g. in Test_2 ladder lane 20) |
| D9 | Well row = **topmost** clear peak of the **median over lanes** | The strongest peak was sometimes the top ladder band |
| D10 | Saturation level = observed maximum if ≥ 5 px reach it | The imager clips at 62846, not 65535 |
| D11 | Size calibration: PCHIP on log10(size) vs. migration; linear extrapolation flagged. Two markers → per-band interpolation across the gel | Monotonic, no overshoot; corrects smiling |
| D12 | Mass curve: volume ∝ ng, through the origin, unsaturated ladder bands only; scaled by ladder loaded / datasheet load | NEB gives ng per band. The user enters only the amount loaded |
| D13 | Ladders in `ladders.json` + user file `~/.autogel/ladders.json`; a `verified` flag per entry | The user wanted a library. Scraping at run time would be fragile, so data is curated and versioned instead |
| D14 | Session JSON stores orientation, params, lanes and bands, plus the log → exact replay | "Log of what was done" must be reproducible, including manual edits |
| D15 | Working image as **float32**, sums in float64 | Exact for 16-bit data, half the RAM of float64 |
| D16 | Okabe-Ito palette as default; band size labels only on marker lanes and the selected lane | Colour-blind users; labels overlapped with 20 lanes (user feedback) |
| D17 | Lane direction (90°) = **weighted sum of four cues**: band shape, periodic well row near the top/bottom, lane-comb periodicity across the gel, edge energy (tie-break). Warn when the sum is near 0 | Band shape alone (which replaced D3's gradient energy in 0.2.1) turned 4 of 95 real gels with heavily smeared lanes by 90°. The well row and the lane comb run across the lanes, so smears don't fool them. Weights checked on the 95 gels both as they are and turned 90° |
| D18 | Ladder library = curated file in the repo + optional downloaded update + ApE + user file; contributions via a GitHub issue form, checked by a maintainer, validated by CI | Users wanted to add ladders and get others' without reinstalling. Scraping vendor pages at run time would be fragile; a reviewed file keeps the data trustworthy |
| D19 | Ladder suggestion from band positions **and** brightness pattern; ladders without a pattern get a neutral brightness score | Positions alone can't tell apart many ladders; without the neutral score, ladders lacking data won |
| D20 | Mass curve: signal per ng as a function of fragment size (from the ladder) | Signal per ng drops 5× from 2 kb to 100 bp on real gels; one slope gave ~90 % typical error, the size curve 19 % |
| D21 | Noise floor = weak relative to the gel **and** near the gel noise | Either criterion alone either removed faint ladder bands (relative) or too little (noise) |
| D22 | Touching bands split by Gaussian areas (on); shoulder splitting optional (off) | Area split fixes biased ratios without changing band counts; shoulder splitting changes counts and misfires on lopsided real bands |
| D23 | Curved lanes estimated per gel (shift + stretch vs. depth), opt-in | Per-lane tracing followed noise and gradients; the gel-wide model is robust on simulations but found a false drift on the lab's straight Test gels |
| D24 | Several channels / added images; default image unchanged; marker lanes may come from another image | Western blots: the chemiluminescence image usually doesn't show the prestained marker; two-colour blots and stain-free images need more than one channel |
| D25 | Lane comparison (band groups, normalisation, RNA) only ADDS columns | Plain densitometry results stay identical for existing users |
| D26 | Learned models only with a fair split (GelGenie train/val vs. test + external labs), as information first, deletions opt-in, never in marker lanes | A model evaluated on its training data looks better than it is; faint ladder bands must never be removed silently |
| D27 | AI second opinion via ONNX + onnxruntime, optional extra, model downloaded on request | PyTorch (~200 MB) would break the "old laptop" goal; nothing leaves the computer |
| D28 | Touching-band peak fit only when R² ≥ 0.97 (was 0.9); background unchanged | On the GelGenie ladder gels the default rolling ball beat every other radius/method; mediocre Gaussian fits of real (tailed) bands cost NEB 1.25 points, Thermo neutral |
| D29 | Marker bands paired with the ladder by local alignment (spacing and brightness per ng change smoothly step to step), classic search kept where clearly smoother; no ladder-guided band splitting/adding | Real ladder lanes: wells, dye front, merged pairs and misses occur together; a global curve bends around them. Splitting/adding bands gave no net gain |

## 4. Test material

| File | Description |
|---|---|
| `<lab share>\Test_1.tif`, `Test_2.tif` | agarose DNA gels; 16-bit LZW TIFF, 2592×1944. 20 lanes, NEB 1 kb Plus (N0550) in lanes 1 and 20. Gel-tray edges and labels visible. Clip level 62846; 30 / 15 saturated px |
| (real-gel tests) | `tests/test_real_gels.py` runs only when the environment variable `AUTOGEL_TEST_GELS` points to the folder with Test_1/Test_2 (the lab share is not in the public repository) |
| `tests/synth.py` | synthetic gels with known positions and intensities (rotation, inversion, dust, clipping, empty lanes) |

## 5. Known limitations / open points
- ng estimates: typical error ~19 % on the Test gels (size-dependent curve, 0.4.0). Part of
  that is the load difference between the two ladder lanes; a per-lane load correction isn't
  modelled.
- Curved-lane tracing (0.4.0) is opt-in: it can find a false drift on straight gels whose
  lower part holds only the ladder lanes. Individual wiggly lanes aren't modelled (only a
  gel-wide fan-out); their outline can be edited by hand.
- Unresolved doublets (shoulders) are only split with the opt-in option, which occasionally
  splits a lopsided single band.
- JPEG / 8-bit images and scratched gels can produce a few spurious faint bands. Use the
  original TIFF, and delete such bands in the viewer.
- The 90° decision is still wrong on 1 of 95 real gels (CB 2022-04-07, a near-tie); the
  *Lane direction uncertain* warning covers it. Rotate 90° fixes it in one click.
- Smiling gels need a marker in the middle for accurate sizes (2–3 % instead of 7–10 %).
- Lane-count auto-estimate without user input can include tray-edge or label peaks next to the
  gel. That's why the GUI asks.
- Ladders not yet verified: GeneRuler 100 bp Plus / 100 bp (sizes and masses) and Bio-Rad
  Precision Plus (sizes). The NEB P7719 protein standard is missing; its band table wasn't
  readable on the product page.
- Ladder suggestion can't tell apart ladders that differ only by a scale factor; the user
  picks from the top 3.

---

## 6. Change log

### 0.5.0.dev2 – 2026-11-19 – branch `develop`: band checker and AI second opinion
User request: use a (AI) model for difficult bands and gels, both as a small own model and as
an optional add-on. Rejected together with the user: cloud vision models (images would leave
the lab, cost, internet, not reproducible) and generic "segment anything" models (large, slow
without a GPU, not made for gels).

**Data and fairness (D26)**
- GelGenie's own split is used for both models:
  - **training:** `images/` and `val_images/` of the Dunn-lab and stella subsets
  - **test:** `test_images/` of every subset plus all `external_gels`, which come from five
    labs never used for training and were GelGenie's own unseen test set
- A detected band counts as real when a hand-outlined band's centre lies inside it.

**Band checker (`autogel/checker.py`, `tools/train_band_checker.py`)**
- **Model:** logistic regression on 11 band features: height vs. the lane's and the gel's
  noise, height vs. the gel's strong bands and the lane's strongest band, lane coverage, width,
  signal at the same rows in the neighbouring lanes, Rf, bands in the lane, peakedness,
  saturation.
- **Fit:** Newton's method with L2 in numpy, no new dependency.
- **Size:** 1.5 KB (`autogel/models/band_checker.json`, with training and test metrics).
- **Test (unseen):**

| Threshold | Kept | Precision | Real kept |
|---|---|---|---|
| none | 100 % | 48.5 % | 100 % |
| 0.1 | 82 % | 57.9 % | 97.8 % |
| 0.2 | 70 % | 64.1 % | 92.7 % |

  At 0.1:
  - external labs 38 → 48 %
  - lab JPEG 33 → 39 % (all real kept)
  - stella 83 → 96 % precision, but only 88 % of the real bands kept
- **Strongest features:** height vs. the gel's strong bands (+), signal in the neighbouring
  lanes at the same height (−; streaks across the gel), lane coverage (+).
- **Test_1/Test_2:** some faint ladder bands next to bright samples score below 0.1 (3 of 38 on
  Test_1). Hence:
  - always computed as information (`band_probability`, grey display)
  - the filter is opt-in and never applied in marker lanes
  - with the filter at 0.1, Test_1 keeps 19/19 + 19/19 ladder bands and loses 8 sample bands

**AI second opinion (`autogel/ai.py`, optional extra `autogel[ai]`)**
- **Model:** GelGenie Universal U-Net (Dec 2023): Apache-2.0, trained on GelGenie's training
  split only. The "Extended" models saw the whole dataset and are not used.
- **Runtime:** ONNX file (57 MB, SHA-256 verified) run with onnxruntime (CPU, ~17 MB) instead of
  PyTorch (~200 MB). About 1.4 s for a 1.4 MP gel.
- **Preprocessing as GelGenie:** grey, scaled by bit depth, padded to a multiple of 32; images
  are capped at 1536 px on the long side. The mask is carried through autogel's orientation.
- **Integration as a second opinion, the user decides:**
  - unconfirmed bands are flagged (grey, `ai_support`)
  - network regions inside a lane without an autogel band are added (source `ai`, green)
  - the model downloads only after the user agrees; no image leaves the computer
- **Test on 43 unseen gels** (`tools/eval_ai.py`):

| Variant | Found | Precision |
|---|---|---|
| autogel | 0.85 | 0.49 |
| autogel, AI-confirmed bands only | 0.79 | **0.91** |
| autogel + AI-added bands | **0.90** | 0.50 |
| GelGenie network alone | 0.86 | 0.84 |

  - The network is weak on the stella gels (found 0.52), a style it wasn't trained on; GelGenie
    used them only for a separate fine-tuned model.
  - Precision of the network alone may be favoured by its outline style matching the hand
    outlines.
- **Download issue:** curl on Windows failed on HuggingFace's CDN (certificate revocation server
  offline, CRYPT_E_REVOCATION_OFFLINE). autogel downloads with Python and certifi instead.
- **Installers:** `-WithAI` / `--with-ai` options; the installer CI installs with the AI extra.

**Tool comparison against true DNA amounts (`tools/eval_quantitation.py`)**
- **Source:** GelGenie paper (Aquilina et al., Nat Commun 2025, CC-BY): quantitation ladder gels
  where every lane holds a ladder with known band masses. The paper's source data gives
  per-lane errors for manual segmentation, GelAnalyzer and LI-COR Image Studio.
- **Protocol** exactly as in their notebook:
  - per lane, fit mass = a·volume + b on all bands but 5 hidden ones
  - error = |true − predicted| / true, averaged over 80 hidden sets (seed 12)
  - lane 2 of gels 14/16/17 excluded; only their 15 NEB gels
- **Our first run** also included gels 25–28 and 30, which the paper didn't evaluate, and 145
  lanes instead of their 102. That gave a median of 48 %; with their gel list it is 100 lanes.
- **Ladders:** NEB Quick-Load Purple 1 kb Plus (N0550) and Invitrogen 1 Kb Plus (10787018,
  current 100–15,000 bp version, now in the library). The 10787018 per-band masses are not in
  the PDF sheet MAN0000898 (it gives sizes and the 500 ng load); they are shown as an image on
  the product page (20 ng 15,000–2,000 bp, 70 ng at 1,500, 30 ng 1,000–200, 40 ng at 100; per
  0.5 µg) and are the same values the GelGenie code uses.
- **Result** (with the peak-fit threshold R² ≥ 0.97, see below), median (middle half) of the
  per-lane error; both ladders in one run, so the random hidden sets are drawn in a different
  order than in the first NEB-only run (16.1 / 26.7 %):

| Method | NEB lanes | NEB median | Thermo lanes | Thermo median |
|---|---|---|---|---|
| autogel, fully automatic (detection + ladder assignment) | 100 | 22.7 % (26.1 before D29) | 180 | 27.6 % |
| autogel, bands identified via the hand outlines | 33 | 15.1 % | 68 | 10.1 % |
| Manual segmentation (paper) | 102 | 17.8 % | 180 | 24.8 % |
| GelAnalyzer, background corrected (paper) | 102 | 17.5 % | 180 | 27.2 % |
| GelAnalyzer, raw (paper) | 102 | 26.1 % | 180 | 27.1 % |
| Image Studio, 3 background variants (paper) | 102 | 24–32 % | 180 | 26–27 % |

- **Same lanes (fair comparison):** the "identified bands" value only exists for lanes where
  every outlined band maps onto exactly one autogel band, i.e. the cleaner lanes. The paper's
  per-gel sheets list lanes in gel order, so on gels with all lanes listed every method can be
  compared on exactly those lanes (`same_lanes` in the script):

| Method, median on the same lanes | NEB (31 lanes) | Thermo (68 lanes) |
|---|---|---|
| autogel, identified bands | 15.1 % | 10.1 % |
| Manual segmentation | 17.0 % | 13.2 % |
| Manual + rolling ball | 16.6 % | 14.9 % |
| GelAnalyzer, raw | 16.9 % | 11.5 % |
| GelAnalyzer, background corrected | **14.7 %** | **9.5 %** |
| Image Studio, 3 variants | 22–28 % | 12–15 % |

- **Reading:**
  - On the same lanes, autogel's band measurement is better than manual segmentation, raw
    GelAnalyzer and Image Studio, and level with GelAnalyzer with background correction
    (15.1 vs 14.7 %, 10.1 vs 9.5 %). The large gap in the first table (10 vs 25 % on Thermo) is mostly lane
    selection, not a better method.
  - Fully automatic, autogel is level with all published tools on the Thermo gels. On the NEB
    gels it is better than raw GelAnalyzer and Image Studio (22.7 % since D29, 26.1 % before)
    but still trails the manually guided methods (~18 %); the remaining gap is mostly the
    900 bp band swallowed by the bright 1000 bp band and missed faint bands.

**Background and peak fit tuned on the ladder gels (user request: try GelAnalyzer's correction)**
- GelAnalyzer's best variant in the paper is a rolling ball on the lane profile; autogel's
  default is already a rolling ball (flat, radius = image height / 40, here 26 px).
- **Bug found first:** `--set` in `eval_quantitation.py` changed the settings after the lane
  profiles were cached, so every setting gave identical numbers. Settings now go in when the
  analysis is created, typed and checked like in the benchmark. No documented number used it.
  The GUI clears the cache before detection and was not affected.
- **Background sweep** (median %, NEB auto / NEB same lanes / Thermo auto / Thermo same lanes):

| Setting | NEB auto | NEB same | Thermo auto | Thermo same |
|---|---|---|---|---|
| default rolling ball, radius ≈ 26 px | **26.4** | **17.3** | **27.5** | 10.1 |
| radius 10 | 34.9 | 20.1 | 29.5 | 19.7 |
| radius 50 | 33.6 | 18.9 | 29.9 | 9.6 |
| radius 100 / 200 | 44–46 | 22–24 | 31 | 14 |
| lane edges | 38.7 | 22.5 | 30.7 | 22.0 |
| no background | 67.5 | 27.2 | 35.4 | 28.0 |

  The default is already the best; no new background option is needed.
- **The real difference was the peak fit of touching bands (D22).** GelAnalyzer integrates the
  profile; autogel fitted Gaussians whenever R² ≥ 0.9. Real ladder bands are not Gaussian
  (tails, tilt), so mediocre fits redistributed signal wrongly. Paired per lane vs. R² ≥ 0.9:

| Peak fit | NEB same lanes | NEB lanes better / worse | Thermo same lanes | Thermo better / worse |
|---|---|---|---|---|
| off | 15.1 % | 30 / 3 | 10.2 % | 18 / 47 (+0.1 points) |
| **R² ≥ 0.97** | **15.1 %** | **25 / 8** (−1.25 points, p < 0.001) | **10.1 %** | 11 / 24 (+0.03 points) |
| R² ≥ 0.99 | 15.1 % | | 10.1 % | |

- **Decision (D28):** `DECONV_MIN_R2 = 0.97` (was 0.9). Clean doublets (simulated test: still
  fitted) keep the benefit; poorly fitting real bands fall back to the valley split. Detection
  is unchanged (the fit only changes volumes).

**Ladder band assignment for the fully automatic route (user request: improve it, D29)**
- **Diagnosis** on the 280 ladder lanes, against the hand outlines:
  - In lanes with complete outlines, a wrong size assignment caused the large errors (NEB:
    median 52 % when wrong vs 14 % when right; 12 of 44 lanes wrong). Measurement was fine.
  - Causes: bands at the wells and the dye front, often several (24_NEB lane 6: all 19 ladder
    bands found plus 4 extras; the old search allowed at most 3 and then paired top -> bottom,
    so every size was wrong); the bright 1000 bp NEB band swallowing 900 bp (22 lanes); a
    merged pair plus a well band giving the ladder's count, so direct pairing shifted the top.
- **Tried and rejected:**
  - global curve + dynamic programming (iterated): too many skips; real migration is not one
    smooth cubic over 100 bp – 10 kb, so true bands looked like outliers
  - beam search over the classic score: too slow (minutes per gel)
  - ladder-guided refinement (split a band where the ladder expects two, add a band where a
    peak sits at an expected size): median NEB 22.7 -> 25.1 %, Thermo 27.4 -> 26.8 %, paired
    not significant. A faint band split off a bright one can't be measured reliably; reverted.
- **Chosen (`markers._align`):** dynamic programming over the last two pairs. Each step costs
  the change of spacing per log(size) against the previous step (measured on true matches:
  ±21 % typical, 95 % within 52 %) and, with masses, the change of brightness per ng (±14 %);
  each unpaired band costs 6, each unpaired size 5; at most 3 skipped in a row. Bands and sizes
  may be skipped together and anywhere. About 0.1 s per lane.
  The classic search is still run and kept where its curve is clearly smoother (by 10 %, or up
  to 1.5x rougher when it uses more bands); this keeps clean, evenly curved ladders (simulated)
  at least as good as before. `suggest_ladders` keeps the classic search (speed).
- **Settings** chosen on assignment precision/recall, not on the quantitation error: a stronger
  brightness weight lowered the error further (NEB 18 %) only by dropping bands that disagree
  with the masses (14.6 of 19 used) - grading itself on easy bands. The chosen setting uses as
  many bands as before.
- **Result** (all 280 lanes; per band against the hand outlines in the 116 complete lanes):

| | NEB precision | NEB recall | NEB median error | Thermo precision | Thermo recall | Thermo median |
|---|---|---|---|---|---|---|
| before | 0.861 | 0.840 | 26.1 % | 0.978 | 0.973 | 27.6 % |
| D29 | **0.979** | **0.946** | **22.7 %** | 0.987 | 0.982 | 27.6 % |

  Paired per lane: NEB 26 better / 9 worse (mean −5.5 points, p < 0.001), Thermo 13 / 11
  (n.s.); bands used 18.0 -> 17.8 (NEB), 16.7 -> 16.6 (Thermo).
- **Synthetic check over the whole library** (83 ladders with 5+ sizes, curved migration): new
  vs. old right: perfect 83/83 vs 83, one band missing 82 vs 82, dye-front band 83 vs 83,
  merged pair 56 vs 54, position noise 83 vs 83. Slowest call 0.2 s (24 sizes).
- **Tests:** wells + dye front beyond 3 extras; merged pair + well band (both fail on the old
  code).
- **Evaluation now as a user:** `eval_quantitation.py` marks every lane as a ladder lane with
  the known ladder and uses `calibrate()`, instead of calling the pairing directly.


### 0.5.0.dev1 – 2026-11-19 – branch `develop`: real-world test data, crop, robustness
User requests: test the new features on real data; use APIs instead of manual downloads; don't
store data in git; tell why a figure set is 5 GB.

**Test data (not in git; `python tools/get_testdata.py` rebuilds it, `--verify` checks it)**
- **Zenodo sets**, checked against Zenodo's MD5 checksums:
  - GelGenie: two more hand-labelled subsets (stella 26 gels, matthew_gels_2 85 gels), 206 gels
    in total
  - Western blots with the authors' densitometry (18338155)
  - Western-blot figure source data (14291079)
- **`tools/oa_testset.py`:** CC-BY open-access papers with raw gel/blot images, via the PLOS
  search API and the Europe PMC REST API.
  - Supplements are found by their captions; data files with the authors' numbers are kept too.
  - Retracted papers are skipped.
  - Currently 49 papers (24 PLOS + 25 Europe PMC; Western blot, RNA, DNA) and 915 pictures.
  - SHA-256 checksums of all 227 files are in `tools/testdata_checksums.json`.
- **Why not via the institute's subscriptions:** that content may not be redistributed (a
  community benchmark must be rebuildable by anyone), and bulk downloads from publisher websites
  can get the whole institute's IP range blocked. Official TDM APIs only, if ever.
- **Why not in git:** GitHub refuses files > 100 MB, repositories should stay < 1–5 GB, Git LFS
  is free only up to 1 GB, and every clone would carry the data forever.
- **The 4.9 GB figure set (14291079):** 13 figure pages at 1200 dpi, CMYK (4 bytes per pixel),
  uncompressed. That's 10,000 × 12,700 px × 4 B ≈ 500 MB per page; the ZIP holds them in 110 MB.
  They're printer files, not raw blot images; the unpacked copies were deleted.

**Robustness (`tools/smoke_images.py`)**
- 238 real-world pictures, **0 crashes**, median 0.8 s per picture.
- **Bug, fixed:** CMYK images (journal figure files) were read inverted (correlation −0.99 with
  the true grey image). They are now converted via Pillow.
- **Memory:** about 50 MB per megapixel; 100-MP figure scans needed up to 7.5 GB. Images over
  25 MP are now reduced by block averaging (relative values within 5 % of full size in the test).

**206 hand-labelled gels**
- autogel finds 85 % of the bands on the 111 new gels too.
- **Benchmark fix:** its orientation reference weighted the hand outlines by area, so a few big
  round blobs or a merged ladder outline could call an upright gel "sideways" (e.g.
  matthew_gels_2/25). Now outlines are counted and only a clear majority is judged: 9 errors on
  196 judged gels (7 distinct; 10 without a clear reference).
- **Orientation remaining errors:** close calls (warning shown) and gels with a ruler and
  notched tray bars in the picture. The ruler's ticks look like a sideways well row, and digits
  are taller than wide.
  - Tried: deciding orientation only on the cropped gel region → 11 errors instead of 9.
    Rejected.

**Western blots vs. the authors (Zenodo 18338155, photographed X-ray films, 36 photos)**
- Default settings: agreement (Spearman ρ) 0.24 (actin) / 0.39 (calponin) / −0.05 (cytokeratin).
- **Cause:** lane placement. The frame, film edge and pen marks add peaks to the lane profile,
  and the lanes on these films are irregular (gaps 45–90 px), so no equidistant comb fits.
- **Measurement itself is fine:** where lanes sit right, duplicate lanes differ by a median of
  6 % in autogel vs. 15 % in the authors' values.
- **Tried: lanes at the n clearest peaks instead of the comb.** No gain on the films, and lab
  gels lost recall (81 → 77 %). Removed.
- **New: `find_gel_region` / auto-crop** (analysed region = the uniform region around the
  centre).
  - On the films: ρ 0.28 / 0.39 / **0.52**; per-patient ratios 0.36 / 0.43 (were 0.23 / 0.17).
  - On the 206 gels: mixed. nathan +4 points, external +1, stella −5.5, matthew −2.7, ladder
    −1.5; 22 gels better and 19 worse, a few by −70 points (part of the gel cut off).
  - **Decision: opt-in** (Expert, *only the gel*), with an editable *analysed region* rectangle
    in the viewer.
- **Checked:** ruff clean; new tests for CMYK, large images, auto-crop / manual crop (core and
  GUI).

### 0.5.0.dev0 – 2026-11-18 – branch `develop`: Western blots, channels, RNA, more ladders
User goal: an open-source community tool for all gel types. The new features are developed on the
branch `develop`; `main` stays at the released 0.4.1 until user feedback is in, and the two are
merged later. 2D gels (spots instead of lanes) are out of scope for now.

**Channels and added images (D24)**
- `load_channels` keeps every channel: RGB (if the colours differ), multi-page TIFF pages.
  The default image is unchanged (colour mean / page 1), so existing results don't change.
- `GelAnalysis.set_channel(k)` analyses another channel; `-1` = all channels, each scaled to its
  own 5th–99.9th percentile and combined by maximum (detects the bands of every colour).
- `add_image(path)` adds a second picture of the same gel/membrane (same raw size, the gel's
  orientation applied, dark-on-light detected per image), e.g. the white-light photo of a blot
  with the prestained marker. Re-attached when a session is loaded.
- `marker_channel`: the marker lanes are read (and their bands detected) from that channel.
- With several channels every band also gets `volume_<channel>`.

**Comparing lanes (`compare.py`, D25)**
- **Band groups:** greedy 1-D clustering across lanes, by log size (±4 %) when calibrated,
  otherwise by Rf (±0.02); at most one band per lane and group. Export: band groups × lanes
  (volume / normalised / fold change).
- **Western blot normalisation:** the loading-control band per lane (strongest band within
  ±10 % of a given size, optionally from another channel) or the total lane signal; `fold_change`
  vs. a reference lane per band group. Two-colour blots: separate channels for loading control
  and target.
- **RNA integrity:** 28S/18S volume ratio (by size with an RNA ladder, else the two strongest
  bands) and the share of signal below 18S.
  - The smear must be measured against the background beside the lane: with the default
    rolling ball, a simulated degraded lane showed only 5 % below 18S (the smear was taken as
    background); with the lane-edge background it shows 85 % (intact: 1 %).

**Ladders, checked on the vendors' pages on 2026-11-18**
- New unit **nt** for RNA.
- Added:
  - NEB ssRNA Ladder N0362: 9000/7000/5000/3000/2000/1000/500 nt, 3000 at double intensity.
    A search summary had suggested a 6000 band; NEB's page says 7000.
  - NEB Low Range ssRNA Ladder N0364: 1000/500/300/150/80/50 nt, 300 at double intensity.
  - Thermo RiboRuler Low Range SM1831: 1000–100.
  - NEB Lambda PFG Ladder N0341: 21 concatemers, 48.5–1018.5 kb.
  - NativeMark LC0725: 1236–20 kDa, 8 bands.
- **Not added:**
  - RiboRuler High Range SM1821: the vendor page gives only the range, no band sizes.
  - NEB Yeast Chromosome PFG Marker N0345: the product page no longer exists; the ApE entry
    stays, unverified.
- NEB's pages refuse automated fetching, so they were read in the built-in browser.

**Tests (simulated, known truth)**
- **Western blot:** chemiluminescence image without visible marker plus a white-light marker
  image (inverted), target 55 kDa at 0.5–4×, loading control 42 kDa with uneven loading.
  - Marker bands 10/10 from the marker image; sizes within 3 %; loading control found in 9/9
    lanes; fold changes within 6 %.
  - The session replay re-attaches the image and gives identical fold changes.
- **Two-colour blot** (target page 1, loading control page 2) with *all channels* detection:
  fold changes within 6 %.
- **RNA:** intact 1.92 (true 2.0) with 1 % below 18S; degraded 1.11 with 85 %.
- **Band groups without calibration** and the band-group matrix; the size check of added images.
- **GUI:** the full Western-blot workflow incl. session restore, and the RNA checkbox.

### 0.4.1 – 2026-11-18 – macOS installer, license, download instructions
User questions: how to get the ZIP, the license, and how Mac users install.

- **License:** BSD-3-Clause, copyright Erik Eppinger. All dependencies are permissive (BSD, MIT,
  PSF, MIT-CMU) except PyQt6 (GPL-3) and Qt (LGPL-3), which are installed by pip on the user's
  machine, not shipped with autogel.
- **Download instructions:** direct ZIP link, the green Code button, Extract All and the
  nested `autogel-main\autogel-main` folder, the Windows "run this file?" prompt.
- **`install_mac.command`** (double-click on macOS; also runs on Linux without the app icon):
  - finds Python ≥ 3.10, preferring 3.12; Apple's `/usr/bin/python3` stub is only asked when
    the Command Line Tools are installed, otherwise it pops up an install dialog
  - opens python.org if no Python is found
  - venv in `~/Library/Application Support/autogel`
  - creates a minimal `~/Applications/autogel.app`; it's made locally, so Gatekeeper doesn't
    block it
  - writes an uninstaller
  - The downloaded `.command` file itself is blocked the first time (unidentified developer);
    the manual explains *Open Anyway* (macOS 15) / right-click → Open (≤ 14) and the
    Terminal fallback.
- **certifi** is now a dependency: python.org's macOS Python has no root certificates, so the
  in-app *Update library* download would fail there.
- **CI `installer check`:**
  - runs the installers on macOS, Ubuntu and Windows exactly as from the ZIP
  - then the installed CLI, a live library update, the full test suite in the installed
    environment and the uninstaller
  - GUI tests need OpenGL: on Linux a virtual screen (xvfb + Mesa); on macOS the runner's window
    server. Qt's "offscreen" mode has no OpenGL, so napari failed to start (first run).
  - **Result (run 44402066431):** macOS, Ubuntu and Windows all pass. On macOS: app created,
    `Info.plist` valid, library update 13 ladders, 110 core + 24 GUI tests passed, uninstaller
    clean.
- **Checked locally:**
  - WSL Ubuntu (Python 3.14 without `venv`): clear error, exit code 1, no half-built
    environment left behind
  - the Windows ZIP route was already checked in 0.4.0

### 0.4.0 – 2026-11-18 – ladder library sharing, ladder recognition, mass curve, noise floor, overlapping and curved lanes
User request: items 1, 2, 3, 4 and 7 of the open-points list, plus a way for users to add ladders
to the library ("a mask for git") that ships with a release or can be updated from the tool.

**Ladder library (D18)**
- **Ladder manager** (*Ladders…*, both modes) replaces the single "save as my ladder" button
  (Expert-only before, name and unit only, no delete, no reference load):
  - new / edit / delete
  - a built-in ladder is never changed: saving an edit creates "… (my copy)"
  - import / export of `.json` files
  - **Suggest for library** opens a pre-filled GitHub issue form (`.github/ISSUE_TEMPLATE/
    ladder.yml`); nothing is sent without the user submitting it
  - **Update library from GitHub** downloads `autogel/ladders.json` from the main branch into
    `~/.autogel/ladders_update.json`
- **Load order:** built-in → update → ApE → user.
- **Validation:** every entry is validated on save, import and download (unit, sizes strictly
  decreasing, as many masses as sizes, reference load). Invalid entries are skipped and never
  crash the program; an invalid or oversized (> 5 MB) download is rejected and the old
  library kept. Files are written atomically.
- **CI:** a GitHub workflow runs `tests/test_ladders.py` whenever a ladder file changes, so a typo
  can't reach users through *Update library*.
- **CLI:** `autogel ladders --import/--replace/--export/--update`.

**Ladder recognition (item 7, D19)**
- `suggest_ladders` pairs every library ladder with the marker lane's bands (like
  `assign_marker`, ±2 bands) and scores:
  - the smoothness of log(size) vs. migration
  - the band count
  - the brightness pattern (datasheet masses or ApE relative intensities)
- Ladders without a brightness pattern get a neutral score; before that, a 25 bp ladder without
  data beat NEB 1 kb Plus on the Test gels.
- **Simulated marker lanes** (156 cases, all distinct library ladders, one band missing in a third):
  top-1 88 %, top-3 99 %, ~0.1 s.
- **Test_1/Test_2:** NEB 1 kb Plus first in all 4 ladder lanes.
- **Limit:** ladders that differ by a scale factor (100 bp vs. 1 kb steps) look alike by position.
- **GUI:** *Suggest ladder* lists the 3 best.

**Mass curve (item 1, D20)**
- **Finding:** on Test_1/Test_2 the signal per ng depends strongly on fragment size, the same way
  in both marker lanes of both gels: highest at 1.5–3 kb, ~20 % of that at 100 bp, ~60 % at
  10 kb. One slope through the origin therefore mis-estimates small fragments 3–5×.
- **New default (`mass_model = "size"`):** log10(volume/ng) is fitted as a polynomial in log
  size (degree 1–3 by the number of distinct sizes) on the unsaturated ladder bands, and each
  sample band is converted at its own size (clamped to the ladder range).
- **Reported:** leave-one-out *typical error* (± factor for a ladder band left out of the fit).

| Gel | one slope: typical error / R² | size-dependent: typical error / R² |
|---|---|---|
| Test_1 | 96 % / 0.67 | 19 % / 0.93 |
| Test_2 | 83 % / 0.76 | 19 % / 0.95 |

- Part of the remaining 19 % is a real load difference between the two ladder lanes (lane 20
  ≈ 25 % brighter).
- New columns: `mass_size_corrected`. The PDF shows signal per ng vs. size; the batch summary
  shows `mass_typical_error`.

**Noise floor (item 2, D21)**
- **Data:** 10 707 bands detected on the 95 benchmark gels, labelled real/spurious by the hand
  masks.
  - Best single feature: height relative to the gel's 95th-percentile band (= *min. height vs.
    gel*). But Test_1's weakest ladder band is 0.8 % of that.
  - Height over the gel-wide noise doesn't depend on bright samples.
- **Rule (default on):** drop a band only if it is < 7 % of the gel's strong bands **and**
  < 10× the gel noise.
  - Benchmark: −28 % of spurious bands for −1.0 % of real ones.
  - The weakest Test ladder band sits at 16.7× noise: margin 1.7×.
  - A stricter version (`< 4 %` and `< 15×`) removed 30 % but with only 11 % margin, so it
    was rejected.

**Overlapping bands (item 3, D22)**
- **Touching bands (default on):** auto bands whose ranges touch are fitted with one Gaussian
  each over the shared range. The shared signal is split by the fitted areas; the fit is
  accepted at R² ≥ 0.9 (0.97 since D28), otherwise the valley split is kept. User-drawn bands are never
  touched.
  - Simulated 4:1 doublet: valley split 4.63, peak fit 4.01.
  - Column `deconvolved`.
- **Shoulders / unresolved doublets (option, default off):** a band is split when its profile
  has two curvature minima (> 6× the curvature noise) **and** two Gaussians fit much better
  than one (RSS < 40 %), the parts are ≥ 1.2 σ apart, and the smaller part carries ≥ 20 % of the
  larger. Saturated bands are skipped (a clipped flat top has two curvature minima).
  - Simulated doublets 1.5:1 to 4:1 at 9–14 px: all found, ratios within 1 %; broad single
    bands not split.
  - On Test_1/2 real ladder bands are slightly lopsided: with a 10 % share, 2 ladder bands per
    gel were split; at 20 %, 1 of 38 is still split on Test_2.
  - Benchmark: strict recall +0.2 to +1.1 points, precision −2 to −4 points. Therefore off.

**Curved lanes (item 4, D23)**
- `Lane.bend` = [[y, dx], …] (centre offset vs. row). It is followed by:
  - profiles (rows grouped by integer offset)
  - coverage, lane-edge background and band re-attachment
  - the viewer (lanes drawn as polygons, bands placed on the bend, polygon edits parsed back)
  - the annotated PNG, size labels and the benchmark boxes
- Old sessions without `bend` load as straight lanes.
- **Tracing, three attempts:**
  1. **Per-lane centroid per strip.** A window of ±0.65 spacing caught the neighbouring ladder
     lane in strips without own bands, giving bends of up to 500 px. With tracking it still
     followed the vignetting gradient and noise.
  2. **Per-lane box matching after removing the slope.** Worse on noisy/dirty gels: random bends
     of 20–60 px.
  3. **Per gel, all lanes at once (kept).** Each of 8 strips is matched to the full-height column
     profile with a shift and a stretch about the gel centre (fan-out / perspective), smoothed
     quadratically over the strips, and accepted only if consistent and > 15 % of the lane
     width.
     - Simulated straight gels (clean, noisy, vignetted, dusty, JPEG, smiling, stacked, side by
       side): no bends.
     - Curved simulation: edge lanes ±23 px (true ~25).
- **But:** on Test_1/Test_2, whose lanes are straight (checked on the image: the bottom ladder
  bands stay centred in the straight boxes), it still found a false drift of up to 30 px. In
  the lower strips only the two ladder lanes carry bands, and two peaks don't constrain a
  stretch. **Decision: tracing is opt-in** (Expert, *follow curved lanes*).
  - Gallery `agarose_curved_lanes`: ratio error 1.9 % straight, 0.7 % traced.

**Skew estimate fixed (found via the benchmark)**
- Since 0.3.1, *Rotate 90°* re-estimates the skew. On CB 2022-04-07 that produced a +1° skew,
  which improved the profile score by only 3 %. The lane comb then locked onto half the
  spacing (recall 0.92 → 0.07).
- `estimate_skew` now applies an angle only if it sharpens the lane profile by ≥ 15 % over 0°;
  a real 2.5° tilt improves it by 168 %.
- Benchmark: 6 images better, none worse.

**Benchmark (default settings, 95 gels):**

| Subset | Found 0.3.1 → 0.4.0 | Strict | Precision 0.3.1 → 0.4.0 |
|---|---|---|---|
| ladder gels | 0.862 → 0.853 | 0.845 | 0.796 → 0.810 |
| lab gels (JPEG) | 0.810 → 0.810 | 0.746 | 0.272 → **0.446** |
| other labs | 0.830 → 0.827 | 0.785 | 0.327 → **0.377** |

With *split shoulders*: strict 0.847 / 0.749 / 0.796, precision 0.766 / 0.427 / 0.354.

- **Also:** the benchmark has `--set name=value` for detection settings.
- **Checked:** ruff clean; `pytest` 136 passed, incl. real gels, the gallery (18 scenarios) and
  the new GUI tests (ladder manager, suggest, all new switches incl. session restore, curved
  lane editing).

### 0.3.1 – 2026-11-18 – lane direction for smeared gels
User request: improve the 90° orientation for gels with heavily smeared lanes; make sure the user
can intervene.

- **Diagnosis:** all 95 benchmark gels are upright; the 4 failures were upright gels turned by
  90°. Their band-shape vote was 0.23–0.59, because smear streaks form tall blobs; the
  edge-energy fallback then decided wrongly.
- **Simply trusting "upright" doesn't work:** tested on the 95 gels plus the same gels turned
  90° (190 cases), an upright-favouring rule only moved the errors to the turned copies
  (9 wrong vs. 8 before).
- **New rule (D17):** a weighted sum of four cues, each positive for "lanes vertical":
  - band shape (vote − 0.5)
  - 0.5 × (periodic well row near top/bottom − same for the turned image)
  - lane comb: periodicity of the mean profile across x − along y
  - ± 0.05 edge energy as tie-break

  **Result:** 188 of 190 cases right, versus 182 for the old rule. Neighbouring weights give
  4–6 errors. About 0.3 s per image.
- **Uncertainty warning:** when |score| < 0.1, the image gets the warning *Lane direction
  uncertain*. This fires on 8 of 190 cases, including both remaining errors. It is not shown
  for an orientation the user set (Rotate 90°, saved session).
- **GUI:** Rotate 90° / Flip ↔ / Flip ↕ are now visible in Easy mode; before, they were
  Expert-only, so an Easy user couldn't fix a wrong orientation. Skew and invert stay in Expert.
- **Bug fix, found in review:** Flip ↔ / Flip ↕ kept the skew angle, but a mirrored image
  tilts the other way, so a gel tilted by 2.5° was left 5° crooked. Flips now negate the angle.
  Rotate 90° re-measures the skew in the new orientation (rotations keep the tilt, so this only
  guards against a skew measured on a wrongly oriented image).
- **Tests:**
  - Flip reverses the skew; Rotate keeps it
  - a smear-fooled band-shape vote on an upright and on a turned synthetic gel is overruled
  - the uncertainty warning
  - rotate/flip visible in Easy mode
  - break test: tiny, constant, pure-noise and 12-MP images; no crash, and featureless images
    score ≈ 0 (warning)
- **Benchmark (default settings):**

| Subset | 90° wrong (0.3.0) | Found (incl. merged) | Strict | Precision |
|---|---|---|---|---|
| ladder gels | 0 (1) | 0.862 | 0.854 | 0.796 |
| lab gels (JPEG) | 0 (0) | 0.810 | 0.749 | 0.272 |
| other labs | 1 (3) | 0.830 | 0.789 | 0.327 |

- **Checked:** ruff clean; `pytest` 98 passed (incl. real gels and gallery).

### 0.3.0 – 2026-11-18 – gels side by side
User request: support two gels lying side by side in one image (seen in the external benchmark
gels).

- **Model:** `Params.n_columns` / `column_splits` and `Lane.column`. Every gel region (gel row ×
  column) gets its own lane comb and well row.
- **Naming:** gels are named A, B, C … in reading order; `gel_column` was added to the band table.
- **Size calibration** uses only markers of the same gel (`tier`, `column`).
- **Automatic split:** the widest low-signal stripe of the band-feature profile within 15–85 %
  of the width. Tray dividers and gel edges are long vertical lines, already removed by the
  vertical top-hat, so they don't count as signal.
- **GUI:** a "gels side by side" setting. The split layer is renamed *gel splits*; horizontal
  lines are gel rows and vertical lines are columns, and both can be dragged.
- **Gallery:** new scenario `agarose_side_by_side` (two 10-lane gels with a bright divider):
  20/20 lanes, 4 × 13/13 ladder bands, 27/27 sample bands, size error ≤ 0.3 %. The split landed
  exactly in the divider gap.
- **Benchmark changes:**
  - The user's layout (lanes per gel, gels side by side) is now estimated in the *oriented*
    frame. Side by side is only assumed for a gap of at least 5 × the lane spacing with
    balanced band clusters (≥ 3 lanes and ≥ 60 % balance per side); the first rule, gap
    > 3 × spacing, also flagged single gels with empty lanes.
  - Orientation errors are now counted separately: if autogel's 90° decision contradicts the
    hand masks, it is recorded and the benchmark continues with the corrected orientation, as a
    user would click Rotate 90°. **4 of 95 gels** had wrong 90° orientation, mostly gels with
    long smeared lanes whose streaks look like bands.
- **Result (default settings, orientation corrected where needed):**

| Subset | 90° wrong | True bands | Found (incl. merged) | Strict | Precision | Inside a lane |
|---|---|---|---|---|---|---|
| ladder gels | 1 | 5421 | 0.868 | 0.860 | 0.798 | 0.915 |
| lab gels (JPEG) | 0 | 625 | 0.810 | 0.749 | 0.272 | 0.902 |
| other labs | 3 | 1297 | 0.831 | 0.784 | 0.327 | 0.938 |

- On the 4 genuine side-by-side gels that autogel oriented correctly, recall stays at 97–100 %
  and precision improves (e.g. 0.22 → 0.36, 0.09 → 0.18).
- **Checked:** ruff clean; `pytest` 93 passed (incl. real gels, the new GUI split test and the
  gallery).

### 0.2.1 – 2026-11-18 – benchmark on 95 hand-labelled real gels
**Data:** GelGenie dataset (Zenodo 10.5281/zenodo.14641949).
- Downloaded with the user's OK: `quantitation_ladder_gels` (35), `nathan_gels` (37, JPEG) and
  `external_gels` (25, five other labs). 97 images, MD5s verified, kept locally in
  `benchmark_data/` (not in the repository).
- Of the 97 images, 95 have usable masks.

**Benchmark (`tools/benchmark_gelgenie.py`):**
- True bands are the connected components of the hand masks.
- The lane count is derived from the true bands (as a user would type it).
- The mask goes through the same orientation autogel applies.
- Reported per gel: found / merged / missed / spurious, recall, precision, and the share of true
  bands inside a placed lane.

**First run, version 0.2.0:**
- Found 79 % (ladder gels) / 68 % (lab JPEG gels) / 75 % (other labs).
- Precision 79 / 22 / 26 %.

**Causes found:**
1. **Orientation.** The "signal centre of mass" rule turned many upright gels upside down, and
   the gradient-energy rule turned some by 90°; the long straight tray and gel edges dominate it.
   New rules:
   - lane direction = area-weighted vote of band-like blobs (wider than tall); lines longer than
     ¼ of the image are ignored
   - upside down only if a periodic well row is clearly stronger at the bottom
     (`row_periodicity`, shared with the gel-row splitter)
   - otherwise the image is left as it is

   Result: 89 of 95 gels kept upright; the simulated rotated gel is still corrected.
2. **Comb spacing bound.** The spacing search had a minimum of W/(n+1)/1.8, which assumes the
   lanes fill the image. On real photos the gel often fills only part of the frame.
   - Removing the bound let narrow combs win inside the bright half of a gel (0_Thermo: 30 px
     instead of 56).
   - A span-based search window broke Test_1/Test_2 (the tray edges widened the span).
   - **Final solution:** a local normalisation of the lane profile (sliding maximum, floored at
     20 % of the global maximum), times a **coverage** factor. The comb must span the signal
     region, measured as the 10th–90th percentile of the profile mass.
3. **Spurious detections:** spread over between bands (1526), below the lowest band (1084),
   above the top band (1019) and outside the band area (631). They are weak: median peak ~1/4 of
   true bands.
   - A gel-wide floor (`min_gel_rel`, 3 % of the gel's 95th-percentile band) raised precision
     to 47 % (lab) / 31 % (external) for about 1 % recall.
   - But on Test_1/2 it removed faint ladder bands (the weakest ladder band is 0.8 % of the 95th
     percentile next to saturated samples).
   - A median-based reference gave almost no gain.
   - **Decision:** the option stays (Expert, "min. height vs. gel"), **default off**.

**Final, default settings:**

| Subset | True bands | Found (incl. merged) | Strict | Precision | True bands inside a lane |
|---|---|---|---|---|---|
| ladder gels (35) | 5420 | 0.886 | 0.874 | 0.802 | 0.933 |
| lab gels, JPEG (37) | 625 | 0.810 | 0.749 | 0.285 | 0.902 |
| other labs (23) | 1297 | 0.798 | 0.753 | 0.239 | 0.924 |

**How to read these numbers:**
- Precision is a lower bound: the annotators did not outline every faint band or smear.
- The lane count was derived from annotated bands, so for gels with many empty lanes it can
  differ from the real number of wells.
- Side-by-side gels are now supported (0.3.0).

**Checked:** ruff clean; `pytest` 91 passed, including the real-gel tests with
`AUTOGEL_TEST_GELS` set; the gallery unchanged; Test_1/2: 20 lanes, 19/19 ladder bands.

### 0.2.0 – 2026-11-18 – logging, GUI test coverage, gel gallery, stacked gels, batch, installer
User requests:
- every GUI option must be tested
- test runs on gels of different quality (simulated if necessary)
- a log format that can be pasted into a lab book
- a visible log line
- a user-selectable size range per gel %
- what does "volume" mean?
- stacked gels with two combs
- batch processing
- an easier installation

**Logging / lab book**
- `report.py`: `describe()` turns every step into a sentence.
- Export now also writes:
  - `<gel>_report.md`: a Markdown report for ELNs with the image SHA-256, settings, step table
    (recalculations merged), per-lane results, warnings and the file list
  - `<gel>_log.jsonl`: machine-readable, one JSON object per step
  - the txt log now shows sentences
- GUI: a status line (⏳ running / ✔ done / ✖ failed), a wait cursor and a collapsible log panel.
  Fed by `GelAnalysis.listeners`; a listener exception can never break the analysis.
- The CLI no longer crashes on cp1252 Windows consoles (µ, ², ⚠).

**GUI test coverage (`tests/test_gui.py`, 19 tests)**
- Every control is exercised on a synthetic gel and must change something measurable: mode
  switch, log, rotate/flip/invert/skew, lane count (given and asked), names, all band-detection
  parameters, re-detect lane, viewer edits, ladder presets, sizes/masses/loads/sample loads, save
  ladder, bad input, palette, label-all, profile navigation and hover, export, session,
  size range, gel rows, batch.
- **Bugs found by these tests and fixed:**
  - Native crash when a second viewer was opened after closing one (timers and connections of the
    closed widget). Added `AutogelWidget.shutdown()` / `closeEvent`.
  - **A band edit in the viewer could be silently dropped.** A spurious "lanes changed" event
    made the edit handler skip band edits. Lane edits now count only if the geometry really
    changed.

**Simulated gel gallery (`tests/gallery.py`, `tests/test_gallery.py`, `tools/run_gallery.py`)**
- 16 scenarios with ground truth: clean, low signal, uneven light, tilted 2.5°, smiling (with 2
  and with 3 markers), dust and scratches, smears and empty lanes, overexposed, 8-bit JPEG,
  rotated plus upside down, 26 narrow lanes, 6 wide lanes, two combs, SDS-PAGE Coomassie RGB
  and greyscale.
- Scored on lanes, ladder bands, sample bands found/missed/extra, size error and relative
  quantification error.
- **Fixes driven by the gallery:**
  - The simulator conserved peak height instead of mass (my error), which caused a fake 11 %
    ratio error. Fixed.
  - Empty lanes and smears made the spacing estimate lock onto 2× the spacing (7 of 12 lanes).
    New `_fit_comb_n`: with the lane count known, an exhaustive comb fit with exactly n teeth,
    then refinement on each lane's plateau.
  - Wide lanes and RGB Coomassie were split in two. The lane-profile baseline window (W/20) was
    narrower than the lanes; now W/4.
  - One extra or missing marker band shifted every size assignment. `assign_marker` now tries
    every way to drop up to 3 surplus bands or ladder sizes and picks the smoothest log(size)
    curve, cross-checking brightness against datasheet masses. The warning names the ignored or
    missing band. On Test_1, deleting the 700 bp band from marker lane 20 is recovered
    automatically.
  - Smiling gels: two edge markers give 7–10 % size error; a middle marker gives 2–3 %.
    Documented as a recommendation.
- **Tried and reverted:** rejecting "bands" that also appear between lanes (meant for JPEG and
  scratch artefacts) also removed real ladder bands on Test_1 (17 and 15 of 19 left).
- **Result:**
  - 0 missed bands in all 16 scenarios
  - median size error ≤ 0.3 % (smiling: 7.4 % with two markers, 2.0 % with three)
  - relative quantification error ≤ 3 % (overexposed: 15 %, flagged as saturated)
  - known: up to 11 spurious faint bands in the JPEG case and 3 in the scratched one

**Valid size range**
- Automatic (whole ladder), suggestions by gel %, or a custom min/max, with 2 % tolerance at the
  limits so a "10 kb" limit includes NEB's 10002 bp band.
- Only ladder bands inside the range feed the size and mass curves; samples outside are flagged
  (`size_in_range`).
- On Test_1 the mass R² goes from 0.67 (whole ladder) to 0.76 (0.5–10 kb) to 0.87 (1–8 kb).
- "Volume" is relabelled as *integrated intensity (a.u.)* in plots and the report.

**Stacked gels / two combs**
- `Params.n_tiers`, `tier_rows`; `Lane.tier`. Lanes, wells and names (A1…, B1…) are found per
  gel row.
- The split is found from the periodicity of well edges along x: ~0.92 for wells against
  0.4–0.64 for band rows. Weighting by edge energy picked a band row instead, so the score uses
  periodicity alone.
- Sizes use only markers of the same gel row.
- GUI: a **gel rows (combs)** setting and a draggable green split line. Moving it re-detects
  the lanes.
- Gallery: 24/24 lanes, 4 × 13/13 ladder bands, 34/34 sample bands.

**Batch (`batch.py`, `autogel batch`, GUI button)**
- Settings come from the current gel or a session file. Each image gets its own orientation,
  lanes, bands, export and session; a failure is logged and the batch continues.
- Combined outputs: `batch_summary.csv`, `batch_bands.csv`, `batch_autogel.xlsx` and
  `batch_report.md` (problem gels first).

**Installation**
- `install_windows.bat` + `tools/install_windows.ps1`:
  - finds Python ≥ 3.10 or offers `winget`
  - creates a venv in `%LOCALAPPDATA%\autogel` and installs autogel with the GUI
  - creates desktop and Start-menu shortcuts, a launcher and an uninstaller
- No admin rights, git or command line needed.
- The test run hit Windows' 260-character path limit (PyQt6 has very deep file names) when
  installing into a long folder. The installer now refuses install paths longer than 70
  characters; the default path is fine.
- Verified: install into a test folder, then import of autogel 0.1.0 and napari 0.10.0, and the
  CLI runs.

**Other**
- `autogel ladders [filter]` lists the library; `--preset` no longer dumps 80 names into the help.
- ApE import: P-numbers (P1473, P9577, P1598) are Sigma, not NEB. Vendor rule fixed and data
  re-imported.
- Real public gel images: GelGenie's Zenodo dataset (575 gels with hand-drawn band masks, CC
  licences) is a good real-world benchmark. Not downloaded yet (needs the user's OK).

### 0.1.3 – 2026-11-18 – readable size labels, visible band boxes
- **Changed** size labels (user feedback: labels sat on the bands and overlapped):
  - New `autogel/labels.py`, shared by the viewer and the exported PNG.
  - Labels sit *beside* the lane, never on the band: marker lanes in the left half are labelled
    on their left, all other lanes on their right. They live on their own point layers
    (`sizes left`, `sizes right`).
  - Sizes are written as 10 kb, 8.0 kb, 1.5 kb, 900.
  - Zoom-aware decluttering: labels closer than one text line at the current zoom are hidden,
    reference bands (3 / 1 / 0.5 kb) first, then by band volume. They reappear when you zoom
    in.
- **Changed** band rectangles to have a light translucent fill. Faint, short bands (11–14 px,
  e.g. Test_1 lane 7 at row 939) showed as a single line when zoomed out, because their top and
  bottom edges merged; bright bands hid one edge (user question).
- **Changed** to use `viewer.scene.camera` (`viewer.camera` is deprecated in napari 0.10).
- **Known**: the selected sample lane's labels can overlap the neighbouring lane's bands, because
  there is no free space between sample lanes.
- **Checked**: ruff clean; `pytest` 52 passed (new label test); rendered viewer at fit and 3× zoom
  and the exported PNG.

### 0.1.2 – 2026-11-18 – GUI usability, ApE ladder import
- **Added** lane navigation in the profile panel: ◀ / ▶ buttons, a lane box, the mouse wheel
  and the arrow keys (user request).
- **Added** hover highlighting: the band under the mouse in the profile is outlined in the
  image, together with its lane. It uses a separate `highlight` layer, so it never interferes
  with edits (user request).
- **Added** Easy / Expert mode, modelled on the simple/expert switch in PrusaSlicer. Easy hides
  orientation, detection parameters, ladder editing, loads and display options. The choice is
  stored with `QSettings` (user feedback: the panel was crowded).
- **Added** `tools/import_ape_ladders.py` and `autogel/ladders_ape.json`:
  - Download approved by the user: `ApE_win_current.zip` (18 MB) from the ApE download page.
  - `ApE.exe` is a Tcl starpack. Its support files were extracted by decompressing the zlib
    streams inside the file; the program was never executed.
  - 67 ladders imported. Repeated sizes become `relative_intensity`, and doublets ≤ 4 % apart
    are merged.
  - Cross-check: ApE's NEB N3200, N3231 and N3232 agree exactly with the datasheet-verified
    sizes. This is now a test.
  - Curated entries take precedence by catalog number.
- **Changed** the profile figure to constrained layout, which stops tight-layout warnings when
  the dock is small.
- **Checked**:
  - ruff clean; `pytest` 51 passed
  - scripted GUI test: lane stepping, clamping at the last lane, hover → highlight (2 shapes)
    and leave → lane only, mode switch and persistence, edits unaffected

### 0.1.1 – 2026-11-18 – ladder library filled from vendor data
- **Changed** `ladders.json` to format 2. Each entry now has:
  - vendor and catalog number
  - `reference_load` (amount, unit, volume, verbatim text)
  - `reference_bands`, `note`, `source`, `verified` and `checked`
- **Verified** (sizes and band masses):
  - NEB N0550, N3200, N3232 and N3231, read from the product pages in the built-in browser;
    NEB blocks automated fetching
  - Thermo GeneRuler SM1331, SM0311 and SM0331, from the user-guide PDFs
    (MAN0013047, MAN0013004, MAN0013012)
- **Verified** (sizes only): PageRuler 26616, 26619 and 26614.
- **Fixed** NEB sizes: they are now the exact values from the datasheets (10002, 8001, …, 2017,
  1517). The NEB 100 bp ladder had 517 and 500 as separate bands; they are one co-migrating
  doublet.
- **Changed** the default preset to "NEB Quick-Load Purple 1 kb Plus (N0550)".
- **Changed** the mass scaling to take its reference load from the library. The GUI enters
  loads in µg and shows each ladder's verification status, load note and source link.
- **Added** `tests/test_ladders.py`: library consistency, and the band masses must add up to the
  stated load. This catches transcription errors: the totals are 1002 ng, 502 ng and 500 ng.
- **Checked**: ruff clean; `pytest` 49 passed (incl. both example gels).

### 0.1.0 – 2026-11-18 – first working version
- **Added** core library `autogel.analysis`: orientation/skew/inversion detection, lane detection,
  lane profiles, three background methods, band detection with dust rejection, size calibration
  (one or two markers), ladder mass curve with load scaling and sample concentration, quantification
  table, session save/load with exact replay, step log.
- **Added** export (`autogel.export`): CSV, XLSX, profile PDF, annotated PNG, log TXT, session JSON.
- **Added** CLI `autogel run|gui`.
- **Added** napari widget with editable lane and band layers, automatic re-quantification after
  edits, lane-profile dock, results table, lane-count confirmation dialog, colour palettes, user
  ladder saving.
- **Added** ladder library (`autogel/ladders.json`); NEB 1 kb Plus checked against the datasheet
  supplied by the user.
- **Added** tests: 23 synthetic/stress tests and 2 regression tests on Test_1/Test_2.
- **Fixed during development** (found by tests or on the real gels):
  - orientation fooled by tray edges (D3)
  - false lanes from tray edges (D4)
  - lane-grid drift and harmonic spacings (D5)
  - dust bands (D8)
  - wrong well row (D9)
  - every band re-tagged as "user" after one edit (an off-by-one between the drawn rectangle and
    the band limits)
  - `set_marker` logged again on every recompute
  - napari 0.10 needs Qt6 and an existing QApplication (switched to PyQt6)
- **Review pass before the first commit**:
  - `ruff` configured (E, F, W, B, I, UP, SIM) and clean
  - removed dead helpers and a sloppy fallback branch
  - lanes dragged outside the image, of zero height or width, or upside down no longer crash
    (break test added)
  - malformed marker or sample-load input now gives a message instead of an exception
  - export errors (e.g. the xlsx still open in Excel) are reported
- **Memory**:
  - the full CLI run on Test_1 peaked at about 1 GB, caused by matplotlib rendering the
    full-resolution annotated PNG
  - the PNG is now drawn with Pillow on an 8-bit copy
  - peak is now **274 MB**, and the whole run including export takes **about 5 s** (D15)
- **Checked**:
  - `pytest`: 28 passed
  - scripted GUI session: open, detect, delete a band, draw a band, move a lane, switch palette,
    take a screenshot
  - both example gels: 20/20 lanes, 19/19 ladder bands in both markers

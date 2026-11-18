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
    environment (Qt offscreen) and the uninstaller
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
  accepted at R² ≥ 0.9, otherwise the valley split is kept. User-drawn bands are never
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

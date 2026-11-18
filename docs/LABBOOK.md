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

## 4. Test material

| File | Description |
|---|---|
| `<lab share>\Test_1.tif`, `Test_2.tif` | agarose DNA gels; 16-bit LZW TIFF, 2592×1944. 20 lanes, NEB 1 kb Plus (N0550) in lanes 1 and 20. Gel-tray edges and labels visible. Clip level 62846; 30 / 15 saturated px |
| `tests/synth.py` | synthetic gels with known positions and intensities (rotation, inversion, dust, clipping, empty lanes) |

## 5. Known limitations / open points
- The ladder mass curve on Test_1 has R² ≈ 0.67. The smallest fragments (100–300 bp) give less
  signal per ng than the datasheet masses predict, and the bright reference bands lie above the
  line. A non-linear (power) mass curve or excluding the extreme bands could be offered.
- Lanes are straight rectangles. Strongly curved lanes need manual adjustment; there is no
  curved-lane model yet.
- Overlapping bands are split at the valley between them, with no Gaussian deconvolution yet.
- Lane-count auto-estimate without user input can include tray-edge or label peaks next to the
  gel. That's why the GUI asks.
- Built-in ladders other than NEB 1 kb Plus are not yet verified against datasheets.

---

## 6. Change log

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

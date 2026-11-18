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
- Ladders not yet verified: GeneRuler 100 bp Plus / 100 bp (sizes and masses) and Bio-Rad
  Precision Plus (sizes). The NEB P7719 protein standard is missing; its band table wasn't
  readable on the product page.
- The ApE ladders have only sizes and relative brightness. They could be used to match ladder
  bands automatically by brightness pattern; that isn't done yet.

---

## 6. Change log

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

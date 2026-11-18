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
- JPEG / 8-bit images and scratched gels can produce a few spurious faint bands. Use the
  original TIFF, and delete such bands in the viewer.
- Smiling gels need a marker in the middle for accurate sizes (2–3 % instead of 7–10 %).
- Lane-count auto-estimate without user input can include tray-edge or label peaks next to the
  gel. That's why the GUI asks.
- Ladders not yet verified: GeneRuler 100 bp Plus / 100 bp (sizes and masses) and Bio-Rad
  Precision Plus (sizes). The NEB P7719 protein standard is missing; its band table wasn't
  readable on the product page.
- The ApE ladders have only sizes and relative brightness. They could be used to match ladder
  bands automatically by brightness pattern; that isn't done yet.

---

## 6. Change log

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

# autogel – User manual

autogel finds lanes and bands on agarose and SDS-PAGE gel images automatically. You can then
correct anything it got wrong. It exports:

- tables of band intensities, sizes and (optionally) DNA amounts
- an intensity profile for every lane
- a complete log that can replay the analysis exactly

It runs on an ordinary office PC or an old laptop. It needs no GPU, and a typical gel uses a few
hundred MB of RAM.

---

## 1. Installation (once)

### Windows – double-click installer (recommended)
1. Download autogel: on the GitHub page, click **Code → Download ZIP**, then unzip the folder
   anywhere, e.g. in *Documents*. If you got the folder from a colleague or a network drive,
   just copy it.
2. Double-click **`install_windows.bat`** inside the folder.
   - If Python 3.10 or newer is missing, the installer offers to install Python 3.12 with
     Windows' own package manager (*winget*). Alternatively, install it from python.org and
     tick *Add python.exe to PATH*.
   - It then downloads autogel's components (about 300 MB, a few minutes) into
     `%LOCALAPPDATA%\autogel`. No administrator rights are needed.
3. Start autogel from the **autogel** icon on the desktop or in the Start menu.

| | |
|---|---|
| Update | download the new version, run `install_windows.bat` again |
| Uninstall | `%LOCALAPPDATA%\autogel\Uninstall autogel.bat` (your own ladders in `%USERPROFILE%\.autogel` are kept) |
| Proxy / no internet | the installer needs internet access for the download; ask IT for the proxy setting (`pip` honours `HTTPS_PROXY`) |

### Developers (command line)

Requires Python ≥ 3.10.

```bash
git clone https://github.com/erikeppinger/autogel
cd autogel
python -m venv .venv
.venv\Scripts\activate          # Windows   (Linux/macOS: source .venv/bin/activate)
pip install -e .[gui,dev]
```

For the command line only, without the napari viewer, run `pip install -e .` instead. This
needs less memory and no OpenGL.

---

## 2. Image requirements – read this first

| Do | Don't |
|---|---|
| Use the **original 16-bit TIFF** from the imager | Use JPG/PNG exports or screenshots. They are compressed and usually auto-contrasted |
| Set the exposure so the bands you want to quantify are **not saturated** | Quantify bands flagged *saturated*. Their intensity is underestimated |
| Image the whole gel including the wells | Crop through the wells |

autogel warns you about JPEG input, 8-bit images and saturated pixels.

---

## 3. Working in the GUI

Start the GUI:

```bash
autogel gui                     # or: autogel gui path\to\gel.tif
```

The **autogel** panel on the right guides you from top to bottom.

**Status line and log:**
- The blue line at the top of the panel always says what autogel is doing (⏳) or has just done (✔).
- **Log ▸** opens the complete list of steps for this gel.

**Easy / Expert mode** (switch at the top of the panel; your choice is remembered):
- **Easy** shows only the essential steps: open the image, set the lane count, detect lanes and
  bands, choose the ladder and marker lanes, quantify, and export.
- **Expert** adds:
  - orientation and skew controls
  - background and detection settings
  - re-detecting bands in a single lane
  - editing ladder sizes and masses
  - datasheet and sample loads
  - display options

  Automatic orientation and the default settings work in both modes.

### 3.1 Image
- **Open image…** loads the gel. autogel automatically:
  - turns the image so the lanes run vertically with the wells on top
  - corrects a slight tilt (skew)
  - detects whether the bands are bright (EtBr, SYBR) or dark (Coomassie, silver)
- If the orientation is wrong, use **Rotate 90°**, **Flip ↔**, **Flip ↕**, the **dark bands**
  checkbox or the **skew °** box. Each correction is logged.
- Saturated pixels are shown in orange-red in the colour-blind-safe palette (red in the classic one).

### 3.2 Lanes
1. Enter the **number of lanes**, including the marker lanes. If you leave it at *ask*, autogel
   estimates the number and asks you to confirm or correct it before it continues. This guards
   against tray edges or labels being counted as lanes.
2. Click **Detect lanes**. Lanes are drawn as rectangles; marker lanes are drawn in orange.
3. To correct lanes, select the **lanes** layer in the layer list and use napari's tools:
   - *select* (arrow) to move or resize a lane
   - *delete* to remove a selected lane
   - *rectangle* to draw a new lane

   Changes apply about half a second after you stop editing.
4. Optional: enter **names** separated by commas (e.g. `M, ctrl, s1, s2, …`).
5. **Several gels in one picture**: stacked (a big gel with two combs) and/or **side by side**
   on the tray.
   - **Gels side by side:** set **gels side by side** to 2 (up to 4). The split is placed in the
     widest empty stripe between the gels; tray dividers don't count. Drag the green vertical
     line in the layer *gel splits* if it's wrong.
   - **Stacked gels:
** set **gel rows (combs)** to 2 (up to 4). autogel finds the
     second row of wells and draws the split as a green horizontal line (layer *gel splits*). If
     it's wrong, drag it to just above the lower wells.
   - After you move a split line, lanes and bands are detected again.
   - **number of lanes** means lanes *per gel*.
   - Gels are named A, B, C, … in reading order (left to right, then top to bottom), and their
     lanes A1…An, B1…Bn, numbered on through all gels (lane n+1 = B1). Give marker lanes by
     number, e.g. `1, 20, 21, 40` for markers at both edges of two 20-lane gels.
   - Sizes are calibrated **only from markers of the same gel**, so every gel needs at least one
     marker. The band table has the columns `gel_row` and `gel_column`.

### 3.3 Bands
- **Detect bands (all lanes)**: finds band peaks in every lane's background-corrected profile.
  Settings:
  - *background*: *rolling_ball* is the default and robust. *lane_edges* takes the signal next
    to each lane as background. *none* applies no correction.
  - *min. signal/noise*: raise it to ignore faint bands, lower it to pick up faint ones.
  - *min. lane coverage*: a real band spans the lane width, a dust speck doesn't. Lower it for
    smeared or tilted bands.
  - *min. height vs. gel* (Expert, default off): ignores bands weaker than this fraction of the
    gel's strong bands.
    - About 0.03 removes most noise "bands" in empty lanes of noisy or JPEG images.
    - Leave it off when faint ladder bands sit next to very bright samples, or they are removed
      too.
- **Re-detect bands in selected lane** re-runs detection for one lane, e.g. after you moved it.
- To correct bands, select the **bands** layer:
  - Each rectangle is one band. Its **height is the integration range**.
  - Resize a rectangle to change what is summed, delete it to drop the band, or draw a new
    rectangle inside a lane to add a band.
  - Bands you added or changed are shown in pink; detected bands in yellow.
- The **lane profile** panel at the bottom:
  - **Change lane:** use the ◀ / ▶ buttons, the *lane* box, the **mouse wheel** over the plot,
    or the **arrow keys** after clicking into the plot.
  - **Hover over a peak** to outline that band and its lane in the image (pink).
  - Selecting a lane or band in the image switches the plot to that lane.
- The plot shows:
  - raw signal, background and corrected signal
  - shaded integration ranges
  - band numbers and sizes
  - saturated points

### 3.4 Marker / calibration
1. Choose the **ladder**. The default is *NEB Quick-Load Purple 1 kb Plus (N0550)*. Sizes (and
   masses, if known) are filled in and can be edited. The grey line under the list shows:
   - whether the entry was checked against the vendor's data (✔) or not yet (⚠)
   - which load the datasheet masses refer to
   - notes such as co-migrating doublets
   - a link to the source
2. Enter the **marker lane(s)**, e.g. `1, 20`. With two markers, sizes are interpolated across
   the gel, which corrects "smiling".
3. Optional DNA amount estimate (only for ladders with band masses):
   - **ladder loaded (µg)**: how much ladder DNA you loaded. Leave it at *as datasheet* if you
     loaded the reference amount.
     - Volume → µg = µl × stock concentration.
     - Example: N0550 is 0.1 µg/µl, so 5 µl = 0.5 µg and the datasheet masses are halved.
   - **datasheet load (µg)** is filled in from the library: NEB 1 kb Plus = 1 µg, NEB 1 kb and
     100 bp = 0.5 µg, GeneRuler = 0.5 µg.
   - autogel fits band volume against ng for the ladder bands (unsaturated bands only, line
     through the origin). The table then gets `est_mass_ng`, and the fit's R² is shown.
   - **sample loads** (`lane:µl`, e.g. `2:5, 3:5`) adds `conc_ng_per_ul`.
4. **Valid size range** (optional): every gel resolves only a certain size range well.
   - *automatic (whole ladder)* is the default.
   - Choosing a gel type (e.g. *agarose 1.0 %* → 500–10 000 bp, *SDS-PAGE 12 %* → 15–100 kDa)
     fills in a typical range; edit **min / max** to set your own.
   - Only ladder bands inside the range are used for the size curve and the ng curve.
   - Sample bands outside it are marked `size_in_range = False` and listed as a warning.
   - Example (Test_1, NEB 1 kb Plus): the ng curve improves from R² 0.67 (whole ladder) to 0.76
     (0.5–10 kb) and 0.87 (1–8 kb).
   - The gel-type ranges are common guidelines. Adjust them to your buffer and run conditions.
5. Click **Calibrate & quantify**.
6. If the number of detected marker bands differs from the ladder (dust, or a band too faint or
   run off), autogel tries every possible assignment and picks the one that gives the smoothest
   size curve. Where the ladder has masses, it also checks that bright bands match the bright
   reference bands. The warning says which band it ignored or assumed missing: check it, and
   add or delete marker bands if it guessed wrong. Doublets such as NEB's 500/517 bp count as
   **one** band.
7. **Smiling gels** (outer lanes run slower than the middle): with markers only at the two edges,
   sizes in the middle lanes can be off by 5–10 %. Load a marker in a middle lane as well (e.g.
   lanes `1, 10, 20`). Sizes are then interpolated between neighbouring markers, which brought
   the error down to 2–3 % in the simulated test gel.

Built-in ladders, as of 2026-11-18:

| Ladder | Masses | Verified |
|---|---|---|
| NEB Quick-Load Purple 1 kb Plus (N0550), NEB 1 kb Plus (N3200) | ng per 1 µg | ✔ |
| NEB 1 kb (N3232), NEB 100 bp (N3231) | ng per 0.5 µg | ✔ |
| GeneRuler 1 kb Plus (SM1331), 1 kb (SM0311), DNA Ladder Mix (SM0331) | ng per 0.5 µg | ✔ |
| PageRuler Prestained (26616), Plus Prestained (26619), Unstained (26614) | – | ✔ sizes |
| GeneRuler 100 bp Plus (SM0321), 100 bp (SM0241), Bio-Rad Precision Plus Dual Color (1610374) | – | ⚠ count/range only |
| 67 more DNA ladders imported from ApE (NEB, Invitrogen, Promega, Sigma, …) | – | ⚠ sizes only |

The ApE ladders come from the free plasmid editor ApE by M. Wayne Davis. They have no ng
values, so they're for sizing only. Where a ladder is in both lists, the checked entry is used.
For the NEB ladders that are in both, ApE's sizes match the datasheets exactly.

**Your own ladders:** edit the sizes and masses, then click **Save sizes/masses as my ladder…**.
It is stored in `~/.autogel/ladders.json` and appears in the list from then on.

> ⚠ Ladder data comes from vendor datasheets. Entries marked `"verified": false` in
> `autogel/ladders.json` have not yet been checked. Lots can change, so compare with the
> datasheet of your lot.

### 3.5 Display
- **colours**: *colour-blind safe (Okabe-Ito)* is the default. *classic* uses cyan, lime, gold
  and red.
- **label all bands**: by default, sizes are shown only for the marker lanes and the selected
  lane, so labels don't overlap.

### 3.6 Results and export
- The table updates after every change. Saturated bands are shown in red.
- **Export…** writes the following to a folder of your choice:

| File | Content |
|---|---|
| `<gel>_bands.csv` | one row per band: lane, name, position, Rf, size, volume, background, peak height, % of lane, saturation, lane coverage, est. mass/concentration |
| `<gel>_lanes.csv` | lane geometry and totals |
| `<gel>_profiles.csv` | full profile of every lane (raw, background, corrected) |
| `<gel>_profiles.pdf` | one profile plot per lane, plus the calibration curves |
| `<gel>_annotated.png` | the gel with lanes, bands and sizes drawn on it |
| `<gel>_autogel.xlsx` | all tables plus the log in one workbook |
| `<gel>_report.md` | **lab-book report** (see below) |
| `<gel>_log.txt` | one readable line per step, with the raw settings |
| `<gel>_log.jsonl` | the same log, machine-readable (one JSON object per line) |
| `<gel>_session.json` | complete state for exact replay |

**Lab-book report (`<gel>_report.md`).** It's written in Markdown, which electronic lab
notebooks (eLabFTW, LabArchives, Benchling, OneNote, Obsidian …) accept directly. Paste it, or
convert it to Word or PDF with any Markdown tool. It contains:
- the image file with its **SHA-256 checksum**, which proves exactly which image was analysed
- size, bit depth, saturation level and the orientation that was applied
- all settings: lane count, background, detection thresholds, marker, ladder load, size range
- the ng mass curve (slope, R², number of bands, range)
- every step with its time in plain language; repeated automatic recalculations are merged
  into one line
- per lane: name, number of bands, sizes, total intensity, total ng, saturated bands
- all warnings and a list of the exported files

Example lines:

```
| 2026-11-18 15:25:29 | detect_lanes | Detected 20 lanes, width 59 px, well row 87 (lane spacing 73.6 px) |
| 2026-11-18 15:25:29 | set_bands    | Bands edited (edited in viewer): 0 added/changed, 1 removed/changed, 94 total |
```

The JSON Lines log carries the same steps with every parameter, for scripts or databases.

- **Save session… / Open session…** stores or restores everything, including manual edits.
  If the image has moved, autogel asks you where it is.

---

## 4. Batch: many gels with the same settings

**In the GUI:**
1. Analyse one gel completely: lane count, gel rows, ladder, marker lanes, size range, and any
   detection settings.
2. Click **Batch: more gels with these settings…**.
3. Select the images and an output folder.

autogel then processes every image automatically. Each image still gets its own orientation,
lanes and bands.

Output:

| File | Content |
|---|---|
| `batch_report.md` | overview: settings used, one line per gel with status **✔ ok / ⚠ check / ✖ error**, problem gels first |
| `batch_summary.csv` | the same as a table (lanes, bands, saturated bands, mass-curve R², warnings) |
| `batch_bands.csv` | all bands of all gels, with a `gel` column, for pivot tables and statistics |
| `batch_autogel.xlsx` | both tables in one workbook |
| `<gel>_autogel/` | the complete normal export of every gel, including its **session** |

**To correct one gel of a batch:** open its `<gel>_session.json` with **Open session…**, fix the
lanes or bands, and export again. A broken or unreadable image is listed as an error and the
batch carries on.

## 5. Command line (no GUI)

```bash
autogel run gel.tif --lanes 20 --marker-lanes 1 20
autogel run gel.tif --lanes 20 --marker-lanes 1 20 --gel "agarose 1.0 %"     # size range by gel type
autogel run gel.tif --lanes 20 --marker-lanes 1 20 --size-range 500 10000
autogel run --session gel_autogel/gel_session.json      # reproduce a saved analysis exactly
autogel batch D:\gels\2026-11 --template gel_autogel/gel_session.json --out D:\results\2026-11
autogel ladders NEB                                    # list ladders (filter by text)
```

Options of `run`:

| Option | Meaning |
|---|---|
| `--preset` | ladder name, as listed by `autogel ladders` |
| `--sizes "10000,8000,…"` | your own ladder sizes |
| `--load` / `--ref-load` | ladder µg loaded / µg the datasheet masses refer to |
| `--invert` | dark bands on a light background |
| `--bg rolling_ball\|lane_edges\|none` | background method |
| `--out DIR` | output folder |

---

## 6. Column reference (`bands.csv`)

| Column | Meaning |
|---|---|
| `y_px`, `y0_px`, `y1_px` | peak row and integration limits (oriented image) |
| `rf` | migration relative to the well (0 = well, 1 = bottom of the image) |
| `size_bp` / `size_kDa` | interpolated from the marker(s). `size_extrapolated` = outside the calibrated ladder bands; `size_in_range` = inside the valid size range and not extrapolated |
| `volume` | background-corrected **integrated intensity** of the band in arbitrary camera units: the sum of all pixel values in the band rectangle (80 % of the lane width × the band's rows) minus the background. "Volume" is the usual densitometry term. **Use this one** |
| `raw_volume`, `background` | the two parts of `volume` |
| `pct_of_lane` | share of all band volume in this lane |
| `saturated`, `saturated_px` | band contains clipped pixels. The volume is too low |
| `lane_coverage` | fraction of the lane width the band covers. Low values suggest dust |
| `est_mass_ng`, `mass_out_of_range` | from the ladder mass curve. *Out of range* = outside the ladder masses |
| `conc_ng_per_ul` | `est_mass_ng / µl loaded` |
| `source` | `auto` (detected) or `user` (added or edited by you) |

---

## 7. How accurate is it?

**Simulated gels with known content** (`tools/run_gallery.py`, 16 situations):
- every band found, sizes within 0.3 %, relative amounts within 3 %
- exceptions:
  - overexposed bands are under-quantified (and flagged)
  - a smiling gel needs a middle marker
  - JPEG and scratched images give a few extra faint bands

**95 real gels outlined by hand** (GelGenie dataset, three sources, `tools/benchmark_gelgenie.py`):

| Gels | autogel finds | of its detections are outlined bands |
|---|---|---|
| ladder test gels | 87 % of the bands | 80 % |
| routine lab gels (JPEG) | 81 % | 27 % |
| gels from five other labs (incl. 5 images with gels side by side) | 83 % | 33 % |

- The 90° orientation was wrong on 4 of 95 images, mostly gels with heavily smeared lanes. Check
  that the lanes run top to bottom and click **Rotate 90°** if not; the numbers above assume
  that correction.
- The last column is a lower bound, because faint bands and smears were often not outlined.

**In practice:**
- Always look over the detected bands of a new gel type.
- Delete spurious ones and add missed ones; both take a click, and the report records it.
- For noisy JPEG gels, try *min. height vs. gel* = 0.03.

## 8. Troubleshooting

| Problem | Fix |
|---|---|
| Wrong number of lanes | Enter the real number under **number of lanes** and detect again. Delete or draw lanes by hand |
| Tray edge or label counted as a lane | Same as above. Entering the lane count avoids it |
| Dust counted as a band | Raise **min. lane coverage**, or delete the band |
| Faint band missed | Lower **min. signal/noise**, or draw the band |
| Two bands merged | Draw two rectangles with the right limits instead of one |
| Sizes wrong | Check that the marker band count matches the ladder (warning under *Calibrate*) |
| Mass R² low | Normal for the smallest and largest fragments, which stain less per ng and are more diffuse. Use only bands within the ladder's mass range, and don't quantify saturated bands |
| Two gels in one image, lanes look wrong | Set **gels side by side** or **gel rows (combs)** to 2 and check the green split line |
| Lanes run across the screen / bands look like vertical streaks | **Rotate 90°**: heavily smeared samples can fool the automatic orientation |
| Many gels to do | Use **Batch** (section 4) with one finished gel as the template |
| GUI does not start on an old PC | Use the command line, or update the graphics driver (napari needs OpenGL 2.1) |

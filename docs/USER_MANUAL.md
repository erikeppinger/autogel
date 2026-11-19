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
1. **Download** autogel:
   - direct link: <https://github.com/erikeppinger/autogel/archive/refs/heads/main.zip>, or
   - on the project page <https://github.com/erikeppinger/autogel>, click the green
     **<> Code** button, then **Download ZIP** at the bottom of the menu.

   The file `autogel-main.zip` (about 150 KB) lands in your *Downloads* folder. If you got the
   folder from a colleague or a network drive, just copy it and skip step 2.
2. **Unzip:** right-click `autogel-main.zip` → **Extract All…** → **Extract**. A folder
   `autogel-main` opens; inside it is a second folder `autogel-main` with the program files.
   You can move it anywhere, e.g. to *Documents*.
   - Don't double-click the installer while you are still *inside* the ZIP: Windows then
     copies out only that one file, the installer can't find the rest and stops with an error.
3. Open the inner `autogel-main` folder (the one that contains `install_windows.bat`,
   `README.md` and the folder `autogel`) and double-click **`install_windows.bat`**.
   - Windows may warn about a downloaded file: click **Run**, or **More info → Run anyway**.
   - If Python 3.10 or newer is missing, the installer offers to install Python 3.12 with
     Windows' own package manager (*winget*). Alternatively, install it from python.org and
     tick *Add python.exe to PATH*.
   - It then downloads autogel's components (about 300 MB, a few minutes) into
     `%LOCALAPPDATA%\autogel`. No administrator rights are needed.
4. Start autogel from the **autogel** icon on the desktop or in the Start menu.

| | |
|---|---|
| Update | download the new version, run `install_windows.bat` again |
| Uninstall | `%LOCALAPPDATA%\autogel\Uninstall autogel.bat` (your own ladders in `%USERPROFILE%\.autogel` are kept) |
| Proxy / no internet | the installer needs internet access for the download; ask IT for the proxy setting (`pip` honours `HTTPS_PROXY`) |

### macOS – double-click installer
1. **Python 3.10 or newer** is needed. macOS doesn't bring a usable one, so install
   **Python 3.12** from <https://www.python.org/downloads/macos/> (the *macOS 64-bit universal2
   installer*, for Intel and Apple-silicon Macs). If you start the autogel installer without it,
   it opens that page for you.
2. **Download** `autogel-main.zip` as for Windows (direct link
   <https://github.com/erikeppinger/autogel/archive/refs/heads/main.zip>, or **<> Code →
   Download ZIP**). Safari unzips it automatically into *Downloads*; otherwise double-click the
   ZIP.
3. Open the folder `autogel-main` and double-click **`install_mac.command`**.
   - **The first time, macOS blocks it** ("Apple could not verify …"), because it doesn't come
     from the App Store or a registered developer.
     - **macOS 15 (Sequoia) and newer:** click **Done**, open **System Settings → Privacy &
       Security**, scroll down to the message about `install_mac.command` and click **Open
       Anyway**, then confirm.
     - **macOS 14 and older:** right-click (or Control-click) the file → **Open** → **Open**.
   - **If that doesn't work:** open *Terminal* (Spotlight: "Terminal"), type `bash ` (with a
     space), drag `install_mac.command` into the Terminal window and press Return.
4. A Terminal window shows the progress. The installer:
   - creates its own environment in `~/Library/Application Support/autogel` and downloads
     autogel's components there (about 300 MB, a few minutes; no administrator rights)
   - creates the app **autogel** in the *Applications* folder of your home folder
5. Start **autogel** from Launchpad or Spotlight, or drag it to the Dock. The app is created on
   your Mac, so macOS doesn't block it.

| | |
|---|---|
| Update | download the new version, run `install_mac.command` again |
| Uninstall | double-click `~/Library/Application Support/autogel/Uninstall autogel.command` (your own ladders in `~/.autogel` are kept) |

### Linux
Use the developer installation below (`pip install -e .[gui]`). On Debian/Ubuntu install the
package `python3-venv` first. `install_mac.command` also works on Linux (without the app
icon): `bash install_mac.command`.

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
- **Easy** shows only the essential steps: open the image (with **Rotate 90°** and the flip
  buttons), set the lane count, detect lanes and bands, choose the ladder and marker lanes,
  quantify, and export.
- **Expert** adds:
  - skew and dark-band (invert) controls
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
- If the orientation is wrong, use **Rotate 90°**, **Flip ↔** or **Flip ↕** (both modes), or the
  **dark bands** checkbox and the **skew °** box (Expert). Each correction is logged, and saved
  sessions keep it.
- When the lane direction is a close call, an orange warning says *Lane direction uncertain*.
  Check that the lanes run top to bottom; click **Rotate 90°** if not.
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

   **Only the gel (auto-crop)** (Expert, off by default): for photos where the gel or membrane
   doesn't fill the picture (photographed X-ray film, a frame, the tray, labels or pen marks
   around it).
   - autogel then looks for lanes only inside the gel. The analysed region is shown as a grey
     rectangle (layer *analysed region*); drag or resize it to set it by hand. Lanes and bands are
     then detected again.
   - On photographed Western-blot films it clearly improved agreement with the authors'
     measurements. On ordinary gel photos it helps some and hurts others (it can cut off part of
     a gel), which is why it is off by default.
   - Very large images (over 25 megapixels, e.g. print-resolution figure scans) are analysed at
     reduced resolution to save memory; a note says so. Relative values are unaffected.

   **Curved lanes** (lanes that fan out towards the bottom, or a photo taken at an angle):
   - Tick *follow curved lanes* (Expert) and detect lanes again. autogel estimates for the whole
     gel how far the lanes drift sideways with depth, and draws bent lanes as polygons.
   - Profiles, bands and the exported image follow the bend. A lane is only bent when the drift
     is consistent and larger than 15 % of the lane width.
   - It is **off by default**: on straight gels where the lower part holds only the ladder
     lanes (like the lab's Test gels), it can still find a false drift. Check the result.
   - A bent lane's outline points can be moved one by one to correct it; dragging the whole
     outline moves the lane.
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
  - *automatic noise floor* (Expert, **on** by default): drops a band only if it is both weak
    (below 7 % of the gel's strong bands) **and** close to the gel's noise (below 10× the noise).
    - This removes noise "bands" in empty lanes of noisy or JPEG images.
    - Faint ladder bands next to very bright samples stay, because they are far above the noise.
    - On the 95 benchmark gels it removed about a quarter of the false detections for about 1 %
      of the real bands. Your Test gels keep all 19 ladder bands.
  - *min. height vs. gel* (Expert, default off): a stricter manual cut at this fraction of the
    gel's strong bands. Leave it off when faint ladder bands sit next to very bright samples.
  - *split touching bands by peak fit* (Expert, **on**): when two bands touch (no gap between
    their ranges), their shared signal is divided by fitted peak areas instead of being cut at
    the dip. This matters for a faint band next to a strong one: a 4:1 pair was measured as
    4.6:1 before and 4.0:1 now. Such bands have `deconvolved = True` in the table.
  - *split shoulders* (Expert, **off**): a band with a clear shoulder (two bands too close to
    show a dip) is split into two.
    - On simulated doublets it finds every pair down to 1:5 and measures their ratio within 1 %.
    - Real bands are often slightly lopsided, so it occasionally splits a single band (1 of 38
      ladder bands on Test_2). Turn it on for close doublets and check the result.
- **Band checker** (learned from hand-labelled gels, always on as information):
  - Every band gets a `band_probability` (0–1) of being a real band. Doubtful bands
    (below 0.1) are drawn in grey.
  - *remove unlikely bands* (Expert, off by default) removes bands below the chosen value after
    detection, **never in marker lanes**.
  - On gels it never saw (five other labs plus GelGenie's test images), 0.1 raised the share of
    real bands among the detections from 49 % to 58 % and kept 98 % of the real bands.
  - The model is a few kilobytes and adds no measurable time.
- **AI second opinion** (button, both modes; optional extra):
  - GelGenie's band-recognition network (Aquilina et al., *Nature Communications* 2025,
    Apache-2.0) looks at the gel.
  - Bands it does **not** confirm are drawn in grey (`ai_support = False` in the table). Bands it
    sees in a lane where autogel has none are added in **green** (source `ai`). Check and
    correct them as usual; the step is logged.
  - On 43 gels the network never saw, keeping only confirmed bands raised precision from 49 %
    to 91 % (bands found 85 % → 79 %). The added bands raised the bands found to 90 %.
  - It runs on the normal CPU, about 1–2 s per gel.
  - **Installation:** tick it in the installer (`install_windows.bat -WithAI`,
    `install_mac.command --with-ai`) or `pip install autogel[ai]` (onnxruntime, about 17 MB).
  - **Model:** the model (57 MB) is downloaded only when you first click the button and agree.
    It is stored in your `.autogel` folder and checked against its checksum.
  - **Privacy:** no image ever leaves your computer.
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
1. Choose the **ladder**. The default is *NEB Quick-Load Purple 1 kb Plus (N0550)*.
   - **Suggest ladder** (after detecting bands, with the marker lane entered) ranks the library
     ladders by how well they fit that lane's bands: positions and the brightness pattern (e.g.
     the bright 3 kb / 1 kb / 500 bp bands of the 1 kb Plus ladder). Pick one of the three best.
     - On Test_1 and Test_2, NEB 1 kb Plus ranks first in all four ladder lanes.
     - Ladders that differ only by a factor (100 bp vs. 1 kb steps) can look alike: check the
       suggestion against what you loaded.
   - Sizes (and masses, if known) are filled in and can be edited. The grey line under the list
     shows:
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
   - autogel learns from the ladder how much signal one ng gives **at each fragment size** and
     uses that for the samples (*mass curve: size-dependent*, the default).
     - Why: on real gels the signal per ng is highest at 1.5–3 kb and drops to about 20 % at
       100 bp and about 60 % at 10 kb (small fragments diffuse, very long ones stain and
       transfer less). A single straight line therefore under-estimates small bands 3–5×.
     - On Test_1/Test_2 the typical error for a ladder band left out of the fit dropped from
       83–96 % (one straight line) to 19 %.
     - The panel shows this *typical error*. The table gets `est_mass_ng` and
       `mass_size_corrected` (False if the band's size is unknown; then the average response is
       used).
     - *mass curve: one slope through origin* (Expert) is the old method.
   - **sample loads** (`lane:µl`, e.g. `2:5, 3:5`) adds `conc_ng_per_ul`.
4. **Valid size range** (optional): every gel resolves only a certain size range well.
   - *automatic (whole ladder)* is the default.
   - Choosing a gel type (e.g. *agarose 1.0 %* → 500–10 000 bp, *SDS-PAGE 12 %* → 15–100 kDa)
     fills in a typical range; edit **min / max** to set your own.
   - Only ladder bands inside the range are used for the size curve and the ng curve.
   - Sample bands outside it are marked `size_in_range = False` and listed as a warning.
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
| NEB ssRNA Ladder (N0362), Low Range ssRNA Ladder (N0364), RiboRuler Low Range RNA Ladder (SM1831) — unit nt | – | ✔ sizes |
| NEB Lambda PFG Ladder (N0341), pulsed-field gels | – | ✔ sizes |
| NativeMark Unstained Protein Standard (LC0725), native protein gels | – | ✔ sizes |
| GeneRuler 100 bp Plus (SM0321), 100 bp (SM0241), Bio-Rad Precision Plus Dual Color (1610374) | – | ⚠ count/range only |
| 67 more DNA ladders imported from ApE (NEB, Invitrogen, Promega, Sigma, …) | – | ⚠ sizes only |

The ApE ladders come from the free plasmid editor ApE by M. Wayne Davis. They have no ng
values, so they're for sizing only. Where a ladder is in both lists, the checked entry is used.
For the NEB ladders that are in both, ApE's sizes match the datasheets exactly.

**Your own ladders and the shared library — Ladders…** (both modes) opens the ladder manager:

| To | Do this |
|---|---|
| Add a ladder | **New**, fill in name, unit, sizes (any order, e.g. pasted from the datasheet), optional masses and the load they refer to, datasheet link → **Save as my ladder** |
| Change a built-in ladder (e.g. your lot differs) | select it, edit, **Save as my ladder** → saved as "… (my copy)"; the original stays |
| Delete | select one of *your* ladders → **Delete** (built-in ones can't be deleted) |
| Give a ladder to a colleague | **Export…** → a small `.json` file; they use **Import…** |
| Get it into autogel for everybody | **Suggest for library…** opens a pre-filled GitHub issue in your browser (needs a free GitHub account; nothing is sent until you submit it). After checking against the datasheet it is added to the library and ships with the next release |
| Get new library ladders without reinstalling | **Update library from GitHub…** downloads the current checked list. Your own ladders are kept. If the download fails or the file is invalid, the old library stays |
| Use it | **Use this ladder** (unsaved edits are saved first) |

Command line: `autogel ladders --import file.json` (add `--replace` to overwrite),
`autogel ladders --export "name" file.json`, `autogel ladders --update`.

Your ladders are stored in `~/.autogel/ladders.json`; a downloaded library update in
`~/.autogel/ladders_update.json`. Every entry is checked when it's saved, imported or downloaded
(units, sizes top to bottom, as many masses as sizes); a broken entry is skipped, never crashes
the program.

> ⚠ Ladder data comes from vendor datasheets. Entries marked `"verified": false` in
> `autogel/ladders.json` have not yet been checked. Lots can change, so compare with the
> datasheet of your lot.

### 3.5 Several images, Western blots, RNA gels

**Channels and added images** (section *Image*, both modes)
- A colour image is analysed as the mean of red/green/blue, and a multi-page TIFF as page 1, as
  before. Under **channel** you can choose a single colour or page instead, or **all
  channels**: every channel is scaled to its own brightness and combined, so the bands of every
  colour are found (e.g. two-colour fluorescent blots).
- **Add image…** adds another picture of the *same* gel or membrane, taken from the same camera
  position without cropping, e.g.:
  - the **white-light photo** of a Western blot, which shows the prestained marker that the
    chemiluminescence image doesn't
  - a total-protein stain (stain-free, Ponceau) for normalisation
  autogel applies the gel's orientation to it and recognises dark-on-light images itself.
- **marker lanes from** (section *Marker*): read the marker lanes from that image. Their bands
  are detected there; everything else stays on the analysed image.
- With several channels, the band table gets one volume column per channel (`volume_red`,
  `volume_page_2`, `volume_<image name>`, …).

**Comparing lanes** (section *Compare lanes*)
- **Band groups:** the same band in different lanes gets the same `band_group` number (1 = top),
  matched by size (within 4 %) when the gel is calibrated, otherwise by relative migration. The
  export contains a table *band groups × lanes* (`*_compare_*.csv` and `compare_*` sheets).
- **normalise to** (Western blots):
  - *loading-control band*: enter its **size** (e.g. `42` for actin in kDa). In every lane the
    strongest band within ±10 % of that size is the loading control (`is_loading_control`), and
    every band is divided by it (`normalised`).
  - *total lane signal*: divide by the lane's whole signal, e.g. from an added total-protein
    image chosen under *loading control from* (Expert).
  - **reference lane**: `fold_change` = normalised value ÷ the same band group's value in that
    lane (e.g. the untreated control).
  - Two-colour blots (Expert): *loading control from* = the channel with the loading control,
    *bands measured in* = the channel with the target; detect bands on *all channels*.
- **RNA integrity:** per lane the 28S/18S volume ratio and the share of the signal below 18S.
  - With an RNA ladder (unit nt), 28S and 18S are found by size (3.2–6 kb and 1.5–2.3 kb);
    otherwise the two strongest bands are taken.
  - Intact mammalian total RNA typically shows a ratio around 2 and little signal below 18S;
    degraded RNA a lower ratio and a smear towards small sizes.
  - The smear is measured against the background *beside* the lane, because the rolling-ball
    background would remove it. These numbers are a guide, not a RIN.
  - Results: the line under *Calibrate & quantify* and the lanes table (`rna_28s_18s`,
    `rna_below_18s_pct`).

On simulated gels:
- **Western blot** (chemiluminescence image plus a separate marker photo; target amounts 0.5–4×
  with uneven loading): sizes within 3 %, and fold changes within 6 % of the truth.
- **RNA:** the 28S/18S ratio comes out as 1.92 (true 2.0), and the degradation smear gives 85 %
  below 18S versus 1 % for intact RNA.

### 3.6 Display
- **colours**: *colour-blind safe (Okabe-Ito)* is the default. *classic* uses cyan, lime, gold
  and red.
- **label all bands**: by default, sizes are shown only for the marker lanes and the selected
  lane, so labels don't overlap.

### 3.7 Results and export
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
| `deconvolved` | True if the band touched a neighbour and its volume comes from the peak fit (see 3.3) |
| `band_probability` | band checker: probability that the band is real (see 3.3) |
| `ai_support` | after *AI second opinion*: True if the network confirms the band |
| `volume_<channel>` | with several channels/images: the band's volume in each of them |
| `band_group`, `group_size` | the same band across lanes (see 3.5) |
| `is_loading_control`, `lc_signal`, `normalised`, `fold_change` | Western-blot normalisation (see 3.5) |
| `est_mass_ng`, `mass_out_of_range` | from the ladder mass curve. *Out of range* = outside the ladder masses |
| `mass_size_corrected` | True if the size-dependent response at this band's size was used (False: size unknown, average response) |
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
| ladder test gels | 85 % of the bands | 81 % |
| routine lab gels (JPEG) | 81 % | 45 % |
| gels from five other labs (incl. 5 images with gels side by side) | 83 % | 38 % |

- The automatic 90° decision was wrong on 1 of 95 images (version 0.3.0: 4), and that image
  got the *Lane direction uncertain* warning. Check that the lanes run top to bottom and click
  **Rotate 90°** if not; the numbers above assume that correction.
- The last column is a lower bound, because faint bands and smears were often not outlined.

**Compared with other tools on true DNA amounts:**
- **Data:** the GelGenie paper's ladder gels, where every band's amount is known, with the
  paper's own protocol and published results.
- **Median error of the predicted DNA amount per lane** (NEB 1 kb Plus / Invitrogen 1 Kb Plus):
  - autogel fully automatic: 26 % / 28 %, like raw GelAnalyzer (26 % / 27 %) and Image Studio
    (24–32 % / 26–27 %)
  - on the lanes where all ladder bands were found, compared on exactly the same lanes:
    autogel 15 % / 10 %, manual segmentation 17 % / 13 %, GelAnalyzer with background
    correction 15 % / 10 %, Image Studio 22–28 % / 12–15 %
- **Practical consequence:** check the ladder bands. When all ladder bands are found, autogel's
  numbers are as good as careful manual work.

**In practice:**
- Always look over the detected bands of a new gel type.
- Delete spurious ones and add missed ones; both take a click, and the report records it.
- The automatic noise floor (on by default) already removes most noise "bands" of JPEG gels;
  *min. height vs. gel* = 0.03 removes more, but can also remove faint ladder bands.

## 8. Troubleshooting

| Problem | Fix |
|---|---|
| Wrong number of lanes | Enter the real number under **number of lanes** and detect again. Delete or draw lanes by hand |
| Tray edge or label counted as a lane | Same as above. Entering the lane count avoids it |
| Dust counted as a band | Raise **min. lane coverage**, or delete the band |
| Faint band missed | Lower **min. signal/noise**, or draw the band |
| Two bands merged | Draw two rectangles with the right limits instead of one |
| Sizes wrong | Check that the marker band count matches the ladder (warning under *Calibrate*) |
| ng estimate typical error high | Check that the right ladder and load are set, and that the ladder bands aren't saturated. The size-dependent mass curve (default) already corrects for small fragments giving less signal per ng |
| A real band was split in two | *split shoulders* (Expert) is on: switch it off, or delete one half and resize the other |
| Lanes are drawn bent although they are straight | Untick *follow curved lanes* (Expert), or drag the lane back |
| Two gels in one image, lanes look wrong | Set **gels side by side** or **gel rows (combs)** to 2 and check the green split line |
| Lanes run across the screen / bands look like vertical streaks | **Rotate 90°**: heavily smeared samples can fool the automatic orientation |
| Many gels to do | Use **Batch** (section 4) with one finished gel as the template |
| GUI does not start on an old PC | Use the command line, or update the graphics driver (napari needs OpenGL 2.1) |

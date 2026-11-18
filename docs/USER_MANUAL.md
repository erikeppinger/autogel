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

Requires Python ≥ 3.10.

```bash
git clone <repository-url> autogel
cd autogel
python -m venv .venv
.venv\Scripts\activate          # Windows   (Linux/macOS: source .venv/bin/activate)
pip install -e .[gui]
```

To use only the command line, without the napari viewer, run `pip install -e .` instead.
This needs less memory and no OpenGL.

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

### 3.3 Bands
- **Detect bands (all lanes)**: finds band peaks in every lane's background-corrected profile.
  Settings:
  - *background*: *rolling_ball* is the default and robust. *lane_edges* takes the signal next
    to each lane as background. *none* applies no correction.
  - *min. signal/noise*: raise it to ignore faint bands, lower it to pick up faint ones.
  - *min. lane coverage*: a real band spans the lane width, a dust speck doesn't. Lower it for
    smeared or tilted bands.
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
4. Click **Calibrate & quantify**.
5. If the number of detected marker bands differs from the ladder, you get a warning. Add or
   delete marker bands until the counts match. Doublets such as NEB's 500/517 bp count as
   **one** band.

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
| `<gel>_log.txt` | human-readable log of every step and setting |
| `<gel>_session.json` | complete state for exact replay |

- **Save session… / Open session…** stores or restores everything, including manual edits.
  If the image has moved, autogel asks you where it is.

---

## 4. Command line (batch, no GUI)

```bash
autogel run gel.tif --lanes 20 --marker-lanes 1 20
autogel run *.tif --lanes 20 --marker-lanes 1 20 --load 5 --ref-load 5
autogel run --session gel_autogel/gel_session.json      # reproduce a saved analysis exactly
```

Options: `--preset`, `--sizes "10000,8000,…"`, `--invert` (dark bands), `--bg rolling_ball|lane_edges|none`,
`--out DIR`.

---

## 5. Column reference (`bands.csv`)

| Column | Meaning |
|---|---|
| `y_px`, `y0_px`, `y1_px` | peak row and integration limits (oriented image) |
| `rf` | migration relative to the well (0 = well, 1 = bottom of the image) |
| `size_bp` / `size_kDa` | interpolated from the marker(s). `size_extrapolated` = outside the ladder range, treat with care |
| `volume` | background-corrected integrated intensity (**use this one**) |
| `raw_volume`, `background` | the two parts of `volume` |
| `pct_of_lane` | share of all band volume in this lane |
| `saturated`, `saturated_px` | band contains clipped pixels. The volume is too low |
| `lane_coverage` | fraction of the lane width the band covers. Low values suggest dust |
| `est_mass_ng`, `mass_out_of_range` | from the ladder mass curve. *Out of range* = outside the ladder masses |
| `conc_ng_per_ul` | `est_mass_ng / µl loaded` |
| `source` | `auto` (detected) or `user` (added or edited by you) |

---

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| Wrong number of lanes | Enter the real number under **number of lanes** and detect again. Delete or draw lanes by hand |
| Tray edge or label counted as a lane | Same as above. Entering the lane count avoids it |
| Dust counted as a band | Raise **min. lane coverage**, or delete the band |
| Faint band missed | Lower **min. signal/noise**, or draw the band |
| Two bands merged | Draw two rectangles with the right limits instead of one |
| Sizes wrong | Check that the marker band count matches the ladder (warning under *Calibrate*) |
| Mass R² low | Normal for the smallest and largest fragments, which stain less per ng and are more diffuse. Use only bands within the ladder's mass range, and don't quantify saturated bands |
| GUI does not start on an old PC | Use the command line, or update the graphics driver (napari needs OpenGL 2.1) |

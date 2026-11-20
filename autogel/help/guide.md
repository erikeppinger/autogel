# Welcome

autogel finds the lanes and bands in a picture of a gel or blot, tells you how big each band
is and how much is in it, and saves the results as a table and a picture.

You go through the panel on the right from top to bottom:

1. **Image** – open your picture
2. **Lanes** – autogel marks the lanes
3. **Bands** – autogel marks the bands
4. **Marker** – tell autogel which lane holds your ladder
5. **Results** – save the table and the picture

You can correct anything autogel did, directly in the picture. Nothing is changed in your
original file.

**Easy and Expert:** Easy shows only what you need for most gels. Expert shows more settings.
Switch at the top of the panel; autogel remembers your choice.

**Where am I?** The blue line at the top always says what autogel is doing or has just done.
**Log** lists every step for this gel.

# Good pictures

The better the picture, the better the numbers.

- **Use the original file from the imaging device** (usually a .tif file). Pictures saved as
  .jpg or .png, and screenshots, have lost detail and often have their brightness changed.
- **Avoid overexposed bands.** If a band is pure white in the middle, its amount is too low in
  the results. autogel marks such bands as *saturated*. Take a shorter exposure if you want to
  measure them.
- **Include the wells** (the pockets at the top) in the picture.
- Straight, evenly run lanes give the best results. Smiling or crooked gels work too, but check
  the lanes.

# Step 1: Open your picture

- Click **Open image…** and choose your file.
- autogel turns the picture so the lanes run from top to bottom with the wells at the top, and
  straightens it if it is slightly tilted.
- **If it is the wrong way round**, use **Rotate 90°**, **Flip ↔** or **Flip ↕**.
- Dark bands on a light background (for example Coomassie-stained gels) are recognised
  automatically.

Several pictures of the same gel or membrane (for example the marker under white light and the
bands in the dark) can be combined with **Add image…**.

# Step 2: Lanes

- Type the **number of lanes** (count your ladder lanes too), then click **Detect lanes**.
  If you leave it at *ask*, autogel counts them and asks you to check.
- Each lane is drawn as a box with its number at the top. Lanes can be different widths, for
  example a narrow ladder lane next to wide sample lanes.
- **Two gels in one picture** (side by side or one above the other): set **gels side by side**
  or **gel rows** first.

**To correct a lane:** click **lanes** in the layer list on the left, then use the arrow tool to
move or resize a box, the delete key to remove it, or the rectangle tool to draw a new one.
Changes take effect half a second after you stop.

# Step 3: Bands

- Click **Detect bands (all lanes)**.
- Each band gets a box. Grey boxes are bands autogel is unsure about – have a look at them.
- **AI second opinion** lets a trained helper check the bands and suggest ones that were
  missed (shown in green). It runs on your computer; the picture is not sent anywhere. The
  first time, it asks before downloading what it needs.

**To correct a band:** click **bands** in the layer list, then move, resize, delete or draw
boxes as for lanes. Hover over a peak in the lane graph at the bottom to see which band it is.
Use the arrow keys or the mouse wheel over the graph to go from lane to lane.

# Step 4: Marker (ladder)

With a ladder, autogel can tell you the size of every band – and, if the ladder's datasheet
lists the amount per band, roughly how much is in each band.

- Choose your **ladder** from the list. **Suggest ladder** compares the bands in your marker
  lane with the ladders autogel knows and proposes the best matches.
- Enter the **marker lane(s)**, for example `1` or `1, 12`. With a ladder at both sides, sizes
  in the middle are more accurate.
- Click **Calibrate & quantify**.
- autogel matches the ladder bands with the sizes, even if a band is missing or an extra spot
  was detected. It tells you what it assumed – check that the sizes look right.
- For protein ladders, choose the entry for your gel type (for example *Bis-Tris MOPS*): the
  same ladder runs differently in different gels.

**Your ladder is not in the list?** Click **Ladders…** to add it from its datasheet. You can
also share it so it is included for everyone.

# Step 5: Results

- **Export…** saves:
  - a table with every band (lane, size, amount, brightness) that opens in Excel
  - the picture with lanes and bands drawn on it
  - a short report of everything that was done, so you can repeat it later
- **Save session…** keeps your work including all corrections, so you can open it again with
  **Open session…**.
- **Batch** analyses more pictures with the same settings. Every picture gets the same number
  of lanes and the same marker lanes, so put only gels with the same layout in one batch.
  **Cancel** stops after the gel that is running; the gels done so far are kept.

# Western blots and RNA gels

- **Western blots:** under **Compare lanes** (Expert) you can compare each lane with a
  loading control (for example actin or GAPDH) or with the total signal of the lane, and show
  every lane relative to a reference lane.
- **RNA gels:** tick **RNA integrity** to see the ratio of the two large ribosomal RNA bands and
  how much of the RNA is broken down.

# When something looks wrong

| What you see | What to do |
|---|---|
| Lanes run sideways | Click **Rotate 90°** |
| Too many or too few lanes | Type the right **number of lanes** and click **Detect lanes** again, or correct the boxes by hand |
| Lanes have the wrong width | Resize the boxes by hand (layer **lanes**, arrow tool) |
| Spots or dust found as bands | Delete them (layer **bands**), or in Expert mode raise **min. height vs. gel** to ignore weak bands |
| A faint band is missing | Draw it (layer **bands**, rectangle tool), or try **AI second opinion** |
| Sizes look wrong | Check the chosen ladder and marker lanes; read the message after **Calibrate & quantify** |
| A band is marked *saturated* | Its amount is too low – take a shorter exposure if you need it |
| autogel feels slow | Close other programs; very large pictures are reduced automatically |

# More help

- **Quick tour:** Help ▸ Quick tour shows each part of the panel.
- **Full manual:** Help ▸ Full manual (online) has every setting and how accurate autogel is.
- **Something does not work, or an idea?** Help ▸ Report a problem opens a form on the
  autogel website. Nothing is sent without you seeing it first.

# autogel

Automatic lane and band detection and quantification for **agarose and SDS-PAGE gels**, with an
interactive **napari** editor for corrections. Built on numpy/scipy, no GPU needed.

- Automatic orientation, lane and band detection (incl. curved lanes), using the lane count you confirm
- Background correction (rolling ball / lane edges), dust rejection, automatic noise floor,
  saturation flags; touching bands split by peak fitting
- Sizes from one or two marker lanes; ng estimates from ladder band masses with a
  size-dependent response curve
- Ladder library (curated + ApE + your own): ladder suggestion from the band pattern,
  import/export, "suggest for library" via GitHub, in-app library update
- Edit lanes and bands directly on the image; every step is logged and replayable
- Several gels per image (side by side and/or stacked), batch processing with the same settings
- Export: CSV / XLSX tables, lane profiles (CSV + PDF), annotated image, lab-book report
  (Markdown), JSON Lines log, session file

**Windows, no command line:**
1. Download **[autogel-main.zip](https://github.com/erikeppinger/autogel/archive/refs/heads/main.zip)**
   (or on this page: green **<> Code** button → **Download ZIP**). It lands in *Downloads*.
2. Right-click the ZIP → **Extract All…** → **Extract**. (Don't start the installer from inside
   the ZIP without extracting; it needs the whole folder.)
3. Open the extracted folder `autogel-main` (Windows shows it one level deeper:
   `autogel-mainutogel-main`) and double-click **`install_windows.bat`**. If Windows asks
   whether to run it, click **Run** (or **More info → Run anyway**).
4. Wait a few minutes (about 300 MB download). Then start autogel from the new desktop icon
   *autogel*.

```bash
pip install -e .[gui]                                     # developers
autogel gui                                               # interactive
autogel run gel.tif --lanes 20 --marker-lanes 1 20        # one gel, no GUI
autogel batch gels/ --template gel_session.json --out results/   # many gels
```

Documentation:
- [docs/USER_MANUAL.md](docs/USER_MANUAL.md) — how to use it
- [docs/LABBOOK.md](docs/LABBOOK.md) — idea, design decisions, change log

Tests: `pip install -e .[dev] && pytest`

## License

BSD 3-Clause, see [LICENSE](LICENSE). The dependencies keep their own licenses (mostly BSD/MIT;
the Qt binding PyQt6, installed with the GUI, is GPL-3). Ladder sizes come from vendor
datasheets and from ApE by M. Wayne Davis (credited in the library files).

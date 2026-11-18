# autogel

Automatic lane and band detection and quantification for **agarose and SDS-PAGE gels**, with an
interactive **napari** editor for corrections. Built on numpy/scipy, no GPU needed.

- Automatic orientation, lane and band detection, using the lane count you confirm
- Background correction (rolling ball / lane edges), dust rejection, saturation flags
- Sizes from one or two marker lanes; ng estimates from ladder band masses
- Edit lanes and bands directly on the image; every step is logged and replayable
- Stacked gels (two combs), batch processing of many gels with the same settings
- Export: CSV / XLSX tables, lane profiles (CSV + PDF), annotated image, lab-book report
  (Markdown), JSON Lines log, session file

**Windows, no command line:** download the ZIP, unzip, double-click `install_windows.bat`.
A desktop icon *autogel* appears.

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

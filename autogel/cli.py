"""Command line: batch analysis without GUI, or start the napari GUI.

    autogel run gel.tif --marker-lanes 1 20 --load 5 --ref-load 5
    autogel run --session gel_autogel/gel_session.json     # exact replay incl. manual edits
    autogel gui [gel.tif]
"""
from __future__ import annotations

import argparse
import os
import sys

from .markers import DEFAULT_PRESET, PRESETS, RESOLUTION_GUIDE, parse_sizes


def main(argv=None):
    # Windows consoles default to cp1252: never crash on µ, ², ⚠ in messages
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser(prog="autogel")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="analyse image(s) automatically and export")
    r.add_argument("images", nargs="*")
    r.add_argument("--session", help="replay a saved session json instead of auto-detection")
    r.add_argument("--out", help="output directory (default: <image>_autogel)")
    r.add_argument("--marker-lanes", type=int, nargs="*", default=[])
    r.add_argument("--preset", default=DEFAULT_PRESET, metavar="LADDER",
                   help="ladder name (see 'autogel ladders'); default: %(default)s")
    r.add_argument("--sizes", help="custom ladder sizes, comma separated (overrides preset)")
    r.add_argument("--load", type=float, help="ladder DNA loaded in µg (scales the datasheet band masses)")
    r.add_argument("--ref-load", type=float, help="µg the datasheet masses refer to (default: from the ladder library)")
    r.add_argument("--lanes", type=int, help="expected number of lanes")
    r.add_argument("--size-range", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="valid size range for calibration (e.g. 500 10000); default: whole ladder")
    r.add_argument("--gel", choices=list(RESOLUTION_GUIDE),
                   help="use the typical size range of this gel type (overridden by --size-range)")
    r.add_argument("--invert", action="store_true", help="dark bands on light background (Coomassie)")
    r.add_argument("--bg", default="rolling_ball", choices=["rolling_ball", "lane_edges", "none"])
    b = sub.add_parser("batch", help="analyse many gels with the settings of a template session")
    b.add_argument("inputs", nargs="+", help="image files, folders or wildcards")
    b.add_argument("--template", required=True,
                   help="session .json of a gel analysed with the wanted settings (lanes, ladder, ...)")
    b.add_argument("--out", required=True, help="output folder")
    lad = sub.add_parser("ladders", help="list, import, export or update the ladder library")
    lad.add_argument("filter", nargs="?", default="", help="text to filter by, e.g. NEB or kDa")
    lad.add_argument("--import", dest="import_file", metavar="FILE",
                     help="add the ladders of a .json file to your library")
    lad.add_argument("--replace", action="store_true", help="with --import: replace ladders of the same name")
    lad.add_argument("--export", nargs=2, metavar=("NAME", "FILE"), help="write one ladder to a file")
    lad.add_argument("--update", action="store_true",
                     help="download the newest curated library from the project's GitHub page")
    g = sub.add_parser("gui", help="open the napari editor")
    g.add_argument("image", nargs="?")
    a = ap.parse_args(argv)

    if a.cmd == "gui":
        from .napari_widget import launch
        launch(a.image)
        return

    if a.cmd == "ladders":
        from . import markers
        if a.import_file:
            r = markers.import_ladders(a.import_file, overwrite=a.replace)
            print(f"added {len(r['added'])}, replaced {len(r['replaced'])}, skipped {len(r['skipped'])}"
                  + (f" (exists; use --replace): {', '.join(r['skipped'])}" if r["skipped"] else ""))
            for err in r["errors"]:
                print("problem:", err)
            return
        if a.export:
            if a.export[0] not in PRESETS:
                ap.error(f"unknown ladder {a.export[0]!r} - see 'autogel ladders'")
            markers.export_ladders([a.export[0]], a.export[1])
            print(f"wrote {a.export[1]}")
            return
        if a.update:
            n, path = markers.update_library()
            print(f"library updated: {n} curated ladders ({path})")
            return
        for name, e in PRESETS.items():
            line = f"{name}  [{e['unit']}, {len(e['sizes'])} bands"
            line += (", masses" if e.get("masses") else "") + (", verified" if e.get("verified") else "") + "]"
            if a.filter.lower() in line.lower():
                print(line)
        return
    if a.cmd == "run" and a.preset not in PRESETS:
        ap.error(f"unknown ladder {a.preset!r} - see 'autogel ladders'")
    if a.cmd == "batch":
        from .batch import collect_images, run_batch
        images = collect_images(a.inputs)
        if not images:
            ap.error("no images found")

        def progress(i, n, name):
            print(f"[{i + 1}/{n}] {name}", flush=True)
        summary = run_batch(images, a.template, a.out, progress)
        for _, r in summary.iterrows():
            print(f"{r['status']:6s} {r['gel']}  {r.get('warnings') or ''}")
        print("batch report:", os.path.join(a.out, "batch_report.md"))
        return

    from .analysis import GelAnalysis, Params
    from .export import export_all
    if a.session:
        ga = GelAnalysis.load_session(a.session)
        print("exported to", export_all(ga, a.out))
        return
    if not a.images:
        ap.error("give image(s) or --session")
    for img in a.images:
        p = Params(invert=True if a.invert else None, n_lanes=a.lanes, bg_method=a.bg)
        ga = GelAnalysis(img, p)
        ga.detect_lanes()
        ga.detect_bands()
        if a.marker_lanes:
            if a.sizes:
                ga.set_marker(a.marker_lanes, sizes=parse_sizes(a.sizes), load=a.load, ref_load=a.ref_load)
            else:
                ga.set_marker(a.marker_lanes, preset=a.preset, load=a.load, ref_load=a.ref_load)
        if a.size_range:
            ga.set_size_range(*a.size_range, reason="--size-range")
        elif a.gel:
            ga.set_size_range(*RESOLUTION_GUIDE[a.gel][:2], reason=a.gel)
        out = export_all(ga, a.out if len(a.images) == 1 else None)
        for w in ga.gel.warnings:
            print("WARNING:", w, file=sys.stderr)
        print(f"{img}: {len(ga.lanes)} lanes, {len(ga.bands)} bands -> {out}")


if __name__ == "__main__":
    main()

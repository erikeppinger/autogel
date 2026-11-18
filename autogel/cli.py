"""Command line: batch analysis without GUI, or start the napari GUI.

    autogel run gel.tif --marker-lanes 1 20 --load 5 --ref-load 5
    autogel run --session gel_autogel/gel_session.json     # exact replay incl. manual edits
    autogel gui [gel.tif]
"""
from __future__ import annotations

import argparse
import sys

from .markers import DEFAULT_PRESET, PRESETS, parse_sizes


def main(argv=None):
    ap = argparse.ArgumentParser(prog="autogel")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="analyse image(s) automatically and export")
    r.add_argument("images", nargs="*")
    r.add_argument("--session", help="replay a saved session json instead of auto-detection")
    r.add_argument("--out", help="output directory (default: <image>_autogel)")
    r.add_argument("--marker-lanes", type=int, nargs="*", default=[])
    r.add_argument("--preset", default=DEFAULT_PRESET, choices=list(PRESETS))
    r.add_argument("--sizes", help="custom ladder sizes, comma separated (overrides preset)")
    r.add_argument("--load", type=float, help="amount of ladder loaded (same unit as --ref-load)")
    r.add_argument("--ref-load", type=float, help="ladder amount the datasheet masses refer to")
    r.add_argument("--lanes", type=int, help="expected number of lanes")
    r.add_argument("--invert", action="store_true", help="dark bands on light background (Coomassie)")
    r.add_argument("--bg", default="rolling_ball", choices=["rolling_ball", "lane_edges", "none"])
    g = sub.add_parser("gui", help="open the napari editor")
    g.add_argument("image", nargs="?")
    a = ap.parse_args(argv)

    if a.cmd == "gui":
        from .napari_widget import launch
        launch(a.image)
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
        out = export_all(ga, a.out if len(a.images) == 1 else None)
        for w in ga.gel.warnings:
            print("WARNING:", w, file=sys.stderr)
        print(f"{img}: {len(ga.lanes)} lanes, {len(ga.bands)} bands -> {out}")


if __name__ == "__main__":
    main()

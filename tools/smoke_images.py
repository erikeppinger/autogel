"""Robustness run over arbitrary gel/blot images (no ground truth).

    python tools/smoke_images.py OUT.csv FOLDER [FOLDER ...] [--min-px 300]

Every image (tif/png/jpg/jp2/bmp, short side >= --min-px, masks/ folders skipped) is loaded,
oriented, and lanes and bands are detected with the default settings and an automatic lane count.
Records crashes, run time, peak memory, orientation decision and counts - it shows whether autogel
copes with real-world files, not whether the result is right.
"""
from __future__ import annotations

import os
import sys
import time
import traceback
import tracemalloc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402

from autogel.analysis import GelAnalysis, Params  # noqa: E402

EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".jp2", ".bmp")
Image.MAX_IMAGE_PIXELS = None


def images(folders, min_px):
    for root in folders:
        for dp, dns, fs in os.walk(root):
            dns[:] = [d for d in dns if d not in ("masks", "val_masks", "test_masks", "__MACOSX")]
            for f in sorted(fs):
                if f.lower().endswith(EXT):
                    p = os.path.join(dp, f)
                    try:
                        with Image.open(p) as im:
                            w, h = im.size
                    except Exception:
                        w = h = 0
                    if min(w, h) >= min_px:
                        yield p, w, h


def main(argv):
    out, folders, min_px = argv[0], [a for a in argv[1:] if not a.startswith("--")], 300
    if "--min-px" in argv:
        min_px = int(argv[argv.index("--min-px") + 1])
        folders = [f for f in folders if f != str(min_px)]
    rows = []
    for p, w, h in images(folders, min_px):
        r = {"image": p, "w": w, "h": h}
        t = time.time()
        tracemalloc.start()
        try:
            ga = GelAnalysis(p, Params())
            o = ga.gel.orientation
            ga.detect_lanes()
            ga.detect_bands()
            r.update(ok=True, lanes=len(ga.lanes), bands=len(ga.bands), rot90=o.get("rot90_k"),
                     invert=o.get("invert"), angle=o.get("angle"), vertical_score=o.get("vertical_score"),
                     channels=len(ga.gel.channels) or 1,
                     warnings=" | ".join(w_ for w_ in ga.gel.warnings if "saturated" not in w_)[:200])
        except Exception as e:
            r.update(ok=False, error=f"{type(e).__name__}: {e}"[:300], trace=traceback.format_exc(limit=3)[-600:])
        r["seconds"] = round(time.time() - t, 2)
        r["peak_mb"] = round(tracemalloc.get_traced_memory()[1] / 1e6)
        tracemalloc.stop()
        rows.append(r)
        print(f"{'ok ' if r.get('ok') else 'ERR'} {r['seconds']:6.1f}s {r['peak_mb']:5d}MB lanes {r.get('lanes')} "
              f"bands {r.get('bands')}  {p[-60:]}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    ok = df[df.ok == True]  # noqa: E712
    print(f"\n{len(df)} images, {len(ok)} ok, {len(df) - len(ok)} errors; median {df.seconds.median():.1f}s, "
          f"max {df.seconds.max():.1f}s, max peak {df.peak_mb.max()} MB")


if __name__ == "__main__":
    main(sys.argv[1:])

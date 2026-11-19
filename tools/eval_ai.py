"""Evaluate the AI second opinion (autogel/ai.py, GelGenie U-Net) on gels it was not trained on.

    python tools/eval_ai.py DATA_DIR OUT_CSV

Held-out gels: GelGenie's test_images of every subset + all external_gels (GelGenie's own unseen
test set). Per gel, against the hand outlines (true band = outlined region >= 15 px):
  autogel        default detection
  ai_filter      autogel bands the network confirms (support)
  union          autogel + bands the network offers in lanes where autogel has none
  gelgenie       the network alone: each connected band region = one detection
Matching as in benchmark_gelgenie.py: a true band is found when its centre lies in a detection
box (+3 px); precision = detections containing a true band / all detections.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import benchmark_gelgenie as bg  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import ndimage  # noqa: E402

from autogel import ai  # noqa: E402

SUBSETS = ["quantitation_ladder_gels", "nathan_gels", "external_gels", "stella_gels_for_finetuning", "matthew_gels_2"]


def score(boxes, gt):
    """(true bands found, detections that contain a true band, detections)"""
    found = sum(1 for (y, x) in gt if any(x0 - 3 <= x <= x1 + 3 and y0 - 3 <= y <= y1 + 3
                                          for x0, x1, y0, y1 in boxes))
    good = sum(1 for x0, x1, y0, y1 in boxes if any(x0 - 3 <= x <= x1 + 3 and y0 - 3 <= y <= y1 + 3
                                                    for (y, x) in gt))
    return found, good, len(boxes)


def band_box(ga, b):
    lane = ga.lanes[b.lane]
    d = lane.dx(b.y)
    return lane.x0 + d, lane.x1 + d, b.y0, b.y1


def main(data: str, out: str):
    rows = []
    for sub in SUBSETS:
        for img, mp in bg.pairs(data, sub):
            folder = os.path.basename(os.path.dirname(img))
            if sub != "external_gels" and not folder.startswith("test"):
                continue
            r = bg.evaluate(img, mp)
            ga, m = r["_ga"], r["_mask"]
            gt = [(rg.centroid[0], rg.centroid[1]) for rg in bg.components(m)]
            rev = ai.review(ga)
            auto = [band_box(ga, b) for b in ga.bands]
            filt = [bx for bx, s in zip(auto, rev["support"]) if s]
            union = auto + [band_box(ga, b) for b in rev["missing"]]
            am = ai.oriented_mask(ga)
            lab, _ = ndimage.label(am)
            gg = [(sl[1].start, sl[1].stop - 1, sl[0].start, sl[0].stop - 1)
                  for k, sl in enumerate(ndimage.find_objects(lab), 1) if sl is not None and (lab[sl] == k).sum() >= 15]
            row = {"subset": sub, "image": os.path.basename(img), "true": len(gt)}
            for name, boxes in (("autogel", auto), ("ai_filter", filt), ("union", union), ("gelgenie", gg)):
                f, g, n = score(boxes, gt)
                row.update({f"{name}_found": f, f"{name}_good": g, f"{name}_det": n})
            rows.append(row)
            print(f"{sub[:10]:10s} {row['image'][:30]:30s} true {len(gt):4d}  autogel {row['autogel_found']:4d}/"
                  f"{row['autogel_det']:4d}  filter {row['ai_filter_found']:4d}/{row['ai_filter_det']:4d}  "
                  f"gelgenie {row['gelgenie_found']:4d}/{row['gelgenie_det']:4d}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print("\nheld-out gels:", len(df))
    for grp, g in [("all", df)] + list(df.groupby("subset")):
        t = g.true.sum()
        line = f"{grp:28s} true {t:5d} |"
        for name in ("autogel", "ai_filter", "union", "gelgenie"):
            rec = g[f"{name}_found"].sum() / t
            prec = g[f"{name}_good"].sum() / max(1, g[f"{name}_det"].sum())
            line += f" {name} R {rec:.2f} P {prec:.2f} |"
        print(line)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

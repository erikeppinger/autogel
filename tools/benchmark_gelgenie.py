"""Benchmark autogel against hand-labelled real gels (GelGenie dataset, Zenodo 10.5281/zenodo.14641949).

    python tools/benchmark_gelgenie.py DATA_DIR OUT_DIR [subset ...]

DATA_DIR contains the unzipped subsets (e.g. quantitation_ladder_gels/quantitation_ladder_gels/
images + masks). Masks: 0 = background, 1 = band.

Per image:
* true bands   = connected components of the mask (>= 15 px)
* lane count   = estimated from the true bands' x positions (what a user would type in)
* the mask is carried through the same orientation (rotation / flip / deskew) autogel applied
* found        = true band whose centre lies inside a detected band box (+3 px)
* merged       = further true bands inside an already matched box (two bands seen as one)
* missed       = true bands inside no box
* spurious     = detected boxes containing no true band
Writes OUT_DIR/benchmark.csv, OUT_DIR/summary.md and annotated images of the worst cases.
"""
from __future__ import annotations

import glob
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from skimage import measure  # noqa: E402

from autogel.analysis import GelAnalysis, Params  # noqa: E402
from autogel.imageio import apply_orientation, load_raw  # noqa: E402

SUBSETS = ["quantitation_ladder_gels", "nathan_gels", "external_gels"]


def pairs(data: str, subset: str):
    root = os.path.join(data, subset, subset)
    for img in sorted(glob.glob(os.path.join(root, "*images", "*"))):
        stem = os.path.splitext(os.path.basename(img))[0]
        mdir = os.path.join(os.path.dirname(os.path.dirname(img)),
                            os.path.basename(os.path.dirname(img)).replace("images", "masks"))
        masks = glob.glob(os.path.join(mdir, stem + ".*"))
        if masks:
            yield img, masks[0]


def components(mask: np.ndarray, min_area: int = 15):
    lab = measure.label(mask > 0, connectivity=2)
    return [r for r in measure.regionprops(lab) if r.area >= min_area]


def estimate_lanes(regions) -> int:
    """Lane count from band x-centres: cluster, then fill gaps with the typical spacing."""
    if not regions:
        return 0
    xs = np.sort([r.centroid[1] for r in regions])
    widths = np.median([r.bbox[3] - r.bbox[1] for r in regions])
    clusters = [[xs[0]]]
    for x in xs[1:]:
        if x - np.mean(clusters[-1]) > 0.5 * widths:
            clusters.append([x])
        else:
            clusters[-1].append(x)
    c = np.array([np.mean(k) for k in clusters])
    if len(c) < 2:
        return 1
    d = np.diff(c)
    sp = np.percentile(d, 25)
    return int(sum(max(1, round(x / sp)) for x in d) + 1)


def estimate_layout(regions, width: int) -> tuple[int, int]:
    """(lanes per gel, gels side by side) as a user would enter them: two gels next to each
    other show up as a gap between band clusters far wider than the lane spacing."""
    if not regions:
        return 0, 1
    xs = np.sort([r.centroid[1] for r in regions])
    widths = np.median([r.bbox[3] - r.bbox[1] for r in regions])
    clusters = [[xs[0]]]
    for x in xs[1:]:
        if x - np.mean(clusters[-1]) > 0.5 * widths:
            clusters.append([x])
        else:
            clusters[-1].append(x)
    c = np.array([np.mean(k) for k in clusters])
    if len(c) >= 4:
        d = np.diff(c)
        sp = np.percentile(d, 25)
        k = int(np.argmax(d))
        mid = 0.5 * (c[k] + c[k + 1])
        nl, nr = k + 1, len(c) - k - 1
        balanced = min(nl, nr) >= 3 and min(nl, nr) / max(nl, nr) >= 0.6
        if d[k] > 5.0 * sp and balanced and 0.15 * width < mid < 0.85 * width:
            left = [r for r in regions if r.centroid[1] < mid]
            right = [r for r in regions if r.centroid[1] >= mid]
            return max(estimate_lanes(left), estimate_lanes(right)), 2
    return estimate_lanes(regions), 1


def evaluate(img: str, mask_path: str, layout: bool = True) -> dict:
    mask_raw, _, _ = load_raw(mask_path)
    regions0 = components(np.asarray(mask_raw))
    t = time.time()
    # the user describes the gel as autogel shows it, i.e. after orientation
    probe = GelAnalysis(img, Params())
    o0 = dict(probe.gel.orientation)
    # Orientation check against the hand masks: true bands are wider than tall when the lanes
    # run vertically. If autogel's 90-degree decision is wrong, record it and continue with the
    # corrected orientation - the user would click "Rotate 90" once.
    raw_regions = components(np.asarray(mask_raw))
    wide = sum(r.area for r in raw_regions if (r.bbox[3] - r.bbox[1]) > (r.bbox[2] - r.bbox[0]))
    tall = sum(r.area for r in raw_regions if (r.bbox[2] - r.bbox[0]) > (r.bbox[3] - r.bbox[1]))
    need_k_odd = tall > wide            # lanes horizontal in the original image
    orientation_error = (o0.get("rot90_k", 0) % 2 == 1) != need_k_odd
    if orientation_error:
        o0 = {**o0, "rot90_k": (o0.get("rot90_k", 0) + 1) % 4, "angle": None}
    m0 = apply_orientation((np.asarray(mask_raw) > 0).astype(np.float32), {**o0, "invert": False}, order=0) > 0.5
    regions0 = components(m0)
    if layout:
        n_lanes, n_cols = estimate_layout(regions0, m0.shape[1])
    else:
        n_lanes, n_cols = estimate_lanes(regions0), 1
    ga = GelAnalysis(img, Params(n_lanes=n_lanes or None, n_columns=n_cols), orientation=o0)
    ga.detect_lanes()
    ga.detect_bands()
    secs = time.time() - t
    o = dict(ga.gel.orientation)
    if o.get("invert"):
        o["invert"] = False
    m = apply_orientation((np.asarray(mask_raw) > 0).astype(np.float32), o, order=0) > 0.5
    regions = components(m)
    gt = [(r.centroid[0], r.centroid[1]) for r in regions]
    boxes = [(ga.lanes[b.lane].x0, ga.lanes[b.lane].x1, b.y0, b.y1) for b in ga.bands]
    owner = [-1] * len(gt)
    for gi, (y, x) in enumerate(gt):
        best = None
        for bi, (x0, x1, y0, y1) in enumerate(boxes):
            if x0 - 3 <= x <= x1 + 3 and y0 - 3 <= y <= y1 + 3:
                d = abs(y - 0.5 * (y0 + y1))
                if best is None or d < best[0]:
                    best = (d, bi)
        if best:
            owner[gi] = best[1]
    used = {}
    for gi, bi in enumerate(owner):
        if bi >= 0:
            used.setdefault(bi, []).append(gi)
    found = len(used)
    merged = sum(len(v) - 1 for v in used.values())
    missed = sum(1 for o_ in owner if o_ < 0)
    spurious = len(boxes) - len(used)
    in_lane = sum(1 for (y, x) in gt if any(l.x0 - 3 <= x <= l.x1 + 3 for l in ga.lanes))
    return {"true_bands": len(gt), "orientation_error": bool(orientation_error),
            "lanes_estimated": n_lanes, "gels_side_by_side": n_cols,
            "lanes_found": len(ga.lanes),
            "detected": len(boxes), "found": found, "merged": merged, "missed": missed,
            "spurious": spurious,
            "recall": round((found + merged) / len(gt), 3) if gt else None,
            "strict_recall": round(found / len(gt), 3) if gt else None,
            "precision": round(found / len(boxes), 3) if boxes else None,
            "bands_in_a_lane": round(in_lane / len(gt), 3) if gt else None,
            "orientation": f"k{o.get('rot90_k', 0)} inv{int(bool(ga.gel.orientation.get('invert')))} "
                           f"ang{o.get('angle', 0):+.2f}",
            "seconds": round(secs, 1), "_ga": ga, "_mask": m}


def overlay(r: dict, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ga, m = r["_ga"], r["_mask"]
    a = ga.data
    lo, hi = np.percentile(a, [1, 99.5])
    fig, ax = plt.subplots(figsize=(10, 7.5))
    ax.imshow(a, cmap="gray", vmin=lo, vmax=hi)
    for c in measure.find_contours(m.astype(float), 0.5):
        ax.plot(c[:, 1], c[:, 0], color="#D55E00", lw=0.6)
    for lane in ga.lanes:
        ax.add_patch(plt.Rectangle((lane.x0, lane.y0), lane.width, lane.y1 - lane.y0, fill=False,
                                   ec="#56B4E9", lw=0.5))
    for b in ga.bands:
        lane = ga.lanes[b.lane]
        ax.add_patch(plt.Rectangle((lane.x0, b.y0), lane.width, b.y1 - b.y0, fill=False, ec="#F0E442", lw=0.6))
    ax.set_title(f"{r['image']}  recall {r['recall']}  precision {r['precision']}  "
                 f"(orange = hand mask, yellow = autogel)", fontsize=8)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main(data: str, out: str, subsets: list[str]):
    os.makedirs(out, exist_ok=True)
    rows = []
    for sub in subsets:
        for img, mp in pairs(data, sub):
            name = os.path.basename(img)
            try:
                r = evaluate(img, mp)
                r.update({"subset": sub, "image": name})
            except Exception as e:
                r = {"subset": sub, "image": name, "error": f"{type(e).__name__}: {e}",
                     "trace": traceback.format_exc(limit=2)}
            rows.append(r)
            print(f"{sub[:12]:12s} {name[:45]:45s} recall {r.get('recall')} precision {r.get('precision')} "
                  f"lanes {r.get('lanes_found')}/{r.get('lanes_estimated')} {r.get('error', '')}", flush=True)
    worst = sorted([r for r in rows if r.get("recall") is not None],
                   key=lambda r: (r["recall"] or 0) + (r["precision"] or 0))[:8]
    for k, r in enumerate(worst):
        overlay(r, os.path.join(out, f"worst_{k + 1}_{os.path.splitext(r['image'])[0][:40]}.png"))
    df = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in rows])
    df.to_csv(os.path.join(out, "benchmark.csv"), index=False)
    ok = df[df.get("error").isna()] if "error" in df else df
    with open(os.path.join(out, "summary.md"), "w", encoding="utf-8") as f:
        f.write("# autogel on hand-labelled real gels (GelGenie dataset)\n\n")
        f.write("| Subset | Images | Errors | 90° orientation wrong | True bands | Recall (incl. merged) "
                "| Strict recall | Precision | Bands inside a lane | s / image |\n"
                "|---|---|---|---|---|---|---|---|---|---|\n")
        for sub, g in df.groupby("subset"):
            gg = ok[ok.subset == sub]
            tb = gg.true_bands.sum()
            f.write(f"| {sub} | {len(g)} | {len(g) - len(gg)} | {int(gg.orientation_error.sum())} | {int(tb)} "
                    f"| {(gg.found.sum() + gg.merged.sum()) / tb:.3f} | {gg.found.sum() / tb:.3f} "
                    f"| {gg.found.sum() / max(1, gg.detected.sum()):.3f} "
                    f"| {(gg.bands_in_a_lane * gg.true_bands).sum() / tb:.3f} | {gg.seconds.mean():.1f} |\n")
    with open(os.path.join(out, "summary.md"), encoding="utf-8") as f:
        print(f.read())


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:] or SUBSETS)

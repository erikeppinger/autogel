"""Quantitation accuracy against known DNA masses, compared with published tools.

    python tools/eval_quantitation.py DATA_DIR SOURCE_DATA_DIR [OUT_DIR] [--set name=value ...]

Data: GelGenie quantitation_ladder_gels (Zenodo 14641949): every lane holds a commercial ladder
with known band masses. NEB lanes (Quick-Load Purple 1 kb Plus, N0550) use autogel's verified
masses. Protocol exactly as Aquilina et al., Nat Commun 2025 (GelGenie), Fig. 1c / S2 / S4:
per lane, fit mass = a * volume + b on all bands but 5 hidden ones, predict the hidden masses,
error = |true - predicted| / true * 100, averaged over 80 random hidden sets (seed 12); lane 2 of
gels 14/16/17_NEB excluded. The published per-lane errors of manual segmentation, GelAnalyzer
and Image Studio come from the paper's source data (SOURCE_DATA_DIR = .../source_data/plots).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import benchmark_gelgenie as bg  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import linregress  # noqa: E402

from autogel.markers import PRESETS, assign_marker  # noqa: E402

NEB = PRESETS["NEB Quick-Load Purple 1 kb Plus (N0550)"]
EXCLUDE = {("14_NEB", 2), ("16_NEB", 2), ("17_NEB", 2)}
# the 15 NEB gels the paper evaluated (GelGenie notebook, gel_based_order)
PAPER_GELS = {"13_NEB", "15_NEB", "14_NEB", "16_NEB", "17_NEB", "18_NEB", "19_NEB", "20_NEB", "22_NEB",
              "23_NEB", "24_NEB", "21_NEB", "29_NEB", "31_NEB", "33_NEB"}


def hidden_linreg(volume, mass, num_hide=5, num_reps=80, rng=None):
    """The paper's per-lane evaluation (GelGenie notebook 'Segmentation Analysis', cell 1)."""
    rng = rng or np.random.default_rng(12)
    n = len(volume)
    if n < num_hide + 2:
        return np.nan
    seen, errs = set(), []
    tries = 0
    while len(seen) < num_reps and tries < 20 * num_reps:
        tries += 1
        combo = tuple(sorted(rng.choice(n, num_hide, replace=False)))
        if combo in seen:
            continue
        seen.add(combo)
        keep = [i for i in range(n) if i not in combo]
        fit = linregress([volume[i] for i in keep], [mass[i] for i in keep])
        errs += [abs(mass[i] - (fit.slope * volume[i] + fit.intercept)) / mass[i] * 100 for i in combo]
    return float(np.mean(errs))


def autogel_lanes(data: str, extra: dict) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(12)
    for img, mp in bg.pairs(data, "quantitation_ladder_gels"):
        gel = os.path.splitext(os.path.basename(img))[0]
        if gel not in PAPER_GELS:
            continue
        r = bg.evaluate(img, mp)                       # same lane count as the benchmark (from the masks)
        ga = r["_ga"]
        for k, v in extra.items():
            setattr(ga.params, k, v)
        if extra:
            ga.detect_bands()
        vols = ga._volumes()
        lanes = sorted(range(len(ga.lanes)), key=lambda i: ga.lanes[i].x)
        for li, i in enumerate(lanes, start=1):
            if (gel, li) in EXCLUDE:
                continue
            bl = [(b.y, v) for b, v in zip(ga.bands, vols) if b.lane == i]
            if len(bl) < 8:
                continue
            ys, sizes, masses, warn = assign_marker([y for y, _ in bl], NEB["sizes"], NEB["masses"],
                                                    volumes=[v for _, v in bl])
            vmap = dict(bl)
            vol = [vmap[y] for y in ys]
            err = hidden_linreg(vol, masses, rng=rng)
            # measurement only: which band is which taken from the hand outlines (as the paper's
            # tools started from identified bands); lanes whose outlines hold all ladder bands
            lane = ga.lanes[i]
            comps = sorted([c for c in bg.components(r["_mask"]) if lane.contains(c.centroid[1], c.centroid[0])],
                           key=lambda c: c.centroid[0])
            err_known = np.nan
            if len(comps) == len(NEB["sizes"]):
                vk = []
                for c in comps:
                    yc = c.centroid[0]
                    hit = [v for b, v in zip(ga.bands, vols) if b.lane == i and b.y0 - 2 <= yc <= b.y1 + 2]
                    # the band autogel measured there (two outlines in one box: shared -> not usable)
                    vk.append(hit[0] if len(hit) == 1 else np.nan)
                lane_bands = [(k, b) for k, b in enumerate(ga.bands) if b.lane == i]
                owners = [next((k for k, b in lane_bands if b.y0 - 2 <= c.centroid[0] <= b.y1 + 2), -1)
                          for c in comps]
                if np.isfinite(vk).all() and len(set(owners)) == len(owners):
                    err_known = hidden_linreg(vk, NEB["masses"], rng=rng)
            rows.append({"Gel": gel, "lane": li, "bands": len(ys), "autogel": err, "autogel_known_bands": err_known})
    return pd.DataFrame(rows)


def main(argv):
    extra = {}
    if "--set" in argv:
        for k in [i for i, a in enumerate(argv) if a == "--set"][::-1]:
            key, _, val = argv[k + 1].partition("=")
            extra[key] = {"0": False, "1": True}.get(val, float(val) if val.replace(".", "").isdigit() else val)
            del argv[k:k + 2]
    data, src = argv[0], argv[1]
    out = argv[2] if len(argv) > 2 else None
    ours = autogel_lanes(data, extra)
    pub = pd.read_excel(os.path.join(src, "figure_1c_source_data.xlsx"))
    pub = pub[pub["Ladder"] == "NEB"]
    ims = pd.read_excel(os.path.join(src, "figure_s4_NEB_Image Studio_source_data.xlsx"))
    seg = pd.read_excel(os.path.join(src, "figure_s2_NEB_segmentation_source_data.xlsx"))
    table = {
        "autogel, end to end (detection + assignment)": ours["autogel"].dropna(),
        "autogel, bands identified as in the paper": ours["autogel_known_bands"].dropna(),
        "Manual segmentation (paper)": pub["Segmentation"],
        "Manual seg. + rolling ball (paper)": seg["Manual + Rolling Ball Background Correction"],
        "GelAnalyzer, raw (paper)": pub["GelAnalyzer (Raw)"],
        "GelAnalyzer, background corr. (paper)": pub["GelAnalyzer (Background Corr.)"],
        "Image Studio, raw (paper)": ims["Image Studio (Raw)"],
        "Image Studio, lane background (paper)": ims["Image Studio (Lane Background Corr.)"],
        "Image Studio, median local bg. (paper)": ims["Image Studio (Median Local Background Corr.)"],
    }
    print(f"NEB 1 kb Plus lanes: autogel evaluated {len(ours)} lanes "
          f"(paper: {len(pub)} lanes); extra settings {extra or 'none'}")
    print(f"{'method':42s} {'lanes':>5s} {'median %':>9s} {'mean %':>7s} {'IQR':>13s}")
    for name, s in table.items():
        s = pd.Series(s).dropna()
        q1, q3 = np.percentile(s, [25, 75])
        print(f"{name:42s} {len(s):5d} {s.median():9.1f} {s.mean():7.1f}   {q1:5.1f}-{q3:5.1f}")
    if out:
        os.makedirs(out, exist_ok=True)
        ours.to_csv(os.path.join(out, "autogel_lanes.csv"), index=False)


if __name__ == "__main__":
    main(sys.argv[1:])

"""Quantitation accuracy against known DNA masses, compared with published tools.

    python tools/eval_quantitation.py DATA_DIR SOURCE_DATA_DIR [OUT_DIR] [--set name=value ...]

Data: GelGenie quantitation_ladder_gels (Zenodo 14641949): every lane holds a commercial ladder
with known band masses: NEB Quick-Load Purple 1 kb Plus (N0550) or Invitrogen 1 Kb Plus
(10787018), both from autogel's ladder library.
Protocol exactly as Aquilina et al., Nat Commun 2025 (GelGenie), Fig. 1c / S2 / S4:
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

LADDERS = {  # name in the gel file -> (ladder, name in the paper's source data)
    "NEB": (PRESETS["NEB Quick-Load Purple 1 kb Plus (N0550)"], "NEB"),
    "Thermo": (PRESETS["Invitrogen 1 Kb Plus (10787018, 100-15,000 bp)"], "ThermoFisher"),
}
EXCLUDE = {("14_NEB", 2), ("16_NEB", 2), ("17_NEB", 2)}
# the 30 gels the paper evaluated (GelGenie notebook, gel_based_order; source data of Fig. S2)
PAPER_GELS = {"13_NEB", "15_NEB", "14_NEB", "16_NEB", "17_NEB", "18_NEB", "19_NEB", "20_NEB", "22_NEB",
              "23_NEB", "24_NEB", "21_NEB", "29_NEB", "31_NEB", "33_NEB",
              "0_Thermo", "1_Thermo", "2_Thermo", "3_Thermo", "4_Thermo", "5_Thermo", "6_Thermo", "7_Thermo",
              "8_Thermo", "9_Thermo", "10_Thermo", "11_Thermo", "12_Thermo", "32_Thermo", "34_Thermo"}


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
        lad = LADDERS[gel.split("_", 1)[1]][0]
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
            ys, sizes, masses, warn = assign_marker([y for y, _ in bl], lad["sizes"], lad["masses"],
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
            if len(comps) == len(lad["sizes"]):
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
                    err_known = hidden_linreg(vk, lad["masses"], rng=rng)
            rows.append({"Gel": gel, "ladder": gel.split("_", 1)[1], "lane": li, "bands": len(ys),
                         "autogel": err, "autogel_known_bands": err_known})
    return pd.DataFrame(rows)


def same_lanes(mine, src, paper_name, seg, ims):
    """Fair comparison: every method on the lanes where autogel's 'identified bands' value exists.

    The paper's per-gel sheets list lanes in gel order, so row k of a gel is taken as lane k. That
    is only safe when every gel has all its lanes in the sheet (no excluded lanes), so gels with
    a different row count are skipped.
    """
    ga = pd.read_excel(os.path.join(src, f"figure_s2_{paper_name}_gelanalyzer_source_data.xlsx"))
    sheets = [seg, ga, ims]
    full = set.intersection(*[set(c[c == c.max()].index) for c in (d["Gel"].value_counts() for d in sheets)])
    m = mine.dropna(subset=["autogel_known_bands"])
    m = m[m["Gel"].isin(full)]
    for d in sheets:
        d = d.assign(lane=d.groupby("Gel").cumcount() + 1)
        m = m.merge(d, on=["Gel", "lane"], suffixes=("", "_dup"))
    if m.empty:
        return
    cols = ["autogel_known_bands", "Manual Segmentation", "Manual + Rolling Ball Background Correction",
            "GelAnalyzer (Raw)", "GelAnalyzer (Background Correction)", "Image Studio (Raw)",
            "Image Studio (Lane Background Corr.)", "Image Studio (Median Local Background Corr.)"]
    print(f"  same {len(m)} lanes (where autogel matched every band to a hand outline), median %:")
    for c in cols:
        if c in m:
            print(f"    {c:48s} {m[c].median():5.1f}")


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
    print(f"extra settings {extra or 'none'}")
    for key, (_, paper_name) in LADDERS.items():
        mine = ours[ours["ladder"] == key]
        pub = pd.read_excel(os.path.join(src, "figure_1c_source_data.xlsx"))
        pub = pub[pub["Ladder"] == paper_name]
        ims = pd.read_excel(os.path.join(src, f"figure_s4_{paper_name}_Image Studio_source_data.xlsx"))
        seg = pd.read_excel(os.path.join(src, f"figure_s2_{paper_name}_segmentation_source_data.xlsx"))
        table = {
            "autogel, end to end (detection + assignment)": mine["autogel"].dropna(),
            "autogel, bands identified as in the paper": mine["autogel_known_bands"].dropna(),
            "Manual segmentation (paper)": pub["Segmentation"],
            "Manual seg. + rolling ball (paper)": seg["Manual + Rolling Ball Background Correction"],
            "GelAnalyzer, raw (paper)": pub["GelAnalyzer (Raw)"],
            "GelAnalyzer, background corr. (paper)": pub["GelAnalyzer (Background Corr.)"],
            "Image Studio, raw (paper)": ims["Image Studio (Raw)"],
            "Image Studio, lane background (paper)": ims["Image Studio (Lane Background Corr.)"],
            "Image Studio, median local bg. (paper)": ims["Image Studio (Median Local Background Corr.)"],
        }
        print()
        print(f"{paper_name} lanes: autogel evaluated {len(mine)} lanes (paper: {len(pub)} lanes)")
        print(f"{'method':42s} {'lanes':>5s} {'median %':>9s} {'mean %':>7s} {'IQR':>13s}")
        for name, ser in table.items():
            ser = pd.Series(ser).dropna()
            if ser.empty:
                print(f"{name:42s}     0")
                continue
            q1, q3 = np.percentile(ser, [25, 75])
            print(f"{name:42s} {len(ser):5d} {ser.median():9.1f} {ser.mean():7.1f}   {q1:5.1f}-{q3:5.1f}")
        same_lanes(mine, src, paper_name, seg, ims)
    if out:
        os.makedirs(out, exist_ok=True)
        ours.to_csv(os.path.join(out, "autogel_lanes.csv"), index=False)


if __name__ == "__main__":
    main(sys.argv[1:])

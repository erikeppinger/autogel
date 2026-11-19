"""Compare autogel with the authors' own densitometry on real Western blots (Zenodo 18338155).

    python tools/eval_wb_18338155.py DATA_DIR [OUT_DIR]

Data: photos of X-ray films, one antibody per photo (beta-actin, calponin, pan-cytokeratin),
4 patients x duplicate lanes plus a marker lane; the file name lists the patients in lane order.
The authors' table gives per patient the two duplicate band values (ImageJ) and the
calponin/actin and cytokeratin/actin ratios.

autogel runs with defaults and 9 lanes; the target band is the band group present in most sample
lanes (strongest in total). Reported: correlation of autogel's lane volumes with the authors'
values (per photo and per antibody) and of the per-patient ratios. This measures agreement with
another method (ImageJ by hand), not absolute truth.
"""
from __future__ import annotations

import glob
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from autogel.analysis import GelAnalysis, Params  # noqa: E402

ANTIBODY = {"b-actin": "actin", "calponin": "cnn", "pan-cytokeratin": "krt"}


def authors_table(path: str) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    hdr = raw.index[raw.iloc[:, 1].astype(str).str.strip().str.lower() == "patient"][0]
    df = raw.iloc[hdr + 1:, 1:].copy()
    df.columns = [str(c).strip().lower() for c in raw.iloc[hdr, 1:]]
    df = df[df["patient"].notna()]
    df["patient"] = df["patient"].astype(str).str.split(",").str[0].str.strip().str.upper()
    return df.set_index("patient")


def patients_from_name(name: str) -> list[str]:
    m = re.search(r"proben\s*([0-9A-Za-z ,.]+?),\s*doppel", name, re.I)
    ids = re.split(r"[,.\s]+", m.group(1).strip()) if m else []
    return [i.upper().lstrip("0") or "0" for i in ids if i]


def lane_volumes(path: str) -> tuple[list[float], str]:
    ga = GelAnalysis(path, Params(n_lanes=9, auto_crop=os.environ.get("AUTOGEL_AUTO_CROP") == "1"))
    ga.detect_lanes()
    ga.detect_bands()
    df = ga.quantify()
    lanes = sorted(range(len(ga.lanes)), key=lambda i: ga.lanes[i].x)
    samples = [i + 1 for i in lanes[1:]]                     # lane 1 = marker
    sub = df[df.lane.isin(samples) & df.band_group.notna()]
    if not len(sub):
        return [np.nan] * 8, "no bands"
    g = sub.groupby("band_group").agg(n=("lane", "nunique"), tot=("volume", "sum"))
    best = g.sort_values(["n", "tot"], ascending=False).index[0]
    vol = dict(zip(sub[sub.band_group == best].lane, sub[sub.band_group == best].volume))
    return [vol.get(L, 0.0) for L in samples], f"group {best} in {g.loc[best, 'n']}/8 lanes"


def main(data: str, out: str | None = None):
    tab = authors_table(glob.glob(os.path.join(data, "*.xlsx"))[0])
    tab.index = [i.lstrip("0") or "0" for i in tab.index]
    rows = []
    for f in sorted(glob.glob(os.path.join(data, "**", "*.jpg"), recursive=True)):
        name = os.path.basename(f)
        ab = next((v for k, v in ANTIBODY.items() if name.lower().startswith(k)), None)
        pts = patients_from_name(name)
        if ab is None or len(pts) != 4:
            print("skip", name)
            continue
        vols, note = lane_volumes(f)
        for k, p in enumerate(pts):
            for d in (0, 1):
                a = tab.loc[p, f"{ab}{d + 1}"] if p in tab.index else np.nan
                rows.append({"image": name, "antibody": ab, "patient": p, "dup": d + 1,
                             "autogel": vols[2 * k + d], "authors": float(a) if pd.notna(a) else np.nan,
                             "note": note})
        print(f"{ab:5s} {note:22s} {name[:60]}", flush=True)
    df = pd.DataFrame(rows).dropna(subset=["authors"])
    if out:
        os.makedirs(out, exist_ok=True)
        df.to_csv(os.path.join(out, "lanes.csv"), index=False)
    print("\nAgreement with the authors' ImageJ values (Spearman rho):")
    per_img = df.groupby("image").apply(lambda g: stats.spearmanr(g.autogel, g.authors)[0], include_groups=False)
    for ab, g in df.groupby("antibody"):
        rho = stats.spearmanr(g.autogel, g.authors)[0]
        imgs = per_img[[i for i in per_img.index if i in set(g.image)]]
        print(f"  {ab:6s} all lanes rho {rho:.2f} (n={len(g)}); per photo median rho {imgs.median():.2f}")
    # per patient: mean of duplicates, ratio target / actin (as in the paper)
    pm = df.groupby(["patient", "antibody"])[["autogel", "authors"]].mean().unstack("antibody")
    for t in ("cnn", "krt"):
        if ("autogel", t) in pm and ("autogel", "actin") in pm:
            ra = (pm[("autogel", t)] / pm[("autogel", "actin")]).replace([np.inf, -np.inf], np.nan)
            rb = pm[("authors", t)] / pm[("authors", "actin")]
            ok = ra.notna() & rb.notna()
            print(f"  {t}/actin per patient: rho {stats.spearmanr(ra[ok], rb[ok])[0]:.2f} (n={int(ok.sum())})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)

"""Train and evaluate the band checker (autogel/checker.py) on hand-labelled gels.

    python tools/train_band_checker.py DATA_DIR [--features CSV] [--write-model]

DATA_DIR = benchmark_data/zenodo_14641949 (GelGenie, CC-BY). GelGenie's own split is kept:
  train = images/ + val_images/ of every subset except external_gels
  test  = test_images/ of every subset + all external_gels (five labs never used for training)
Each detected band (default settings) is labelled real if a hand-outlined band's centre lies in
its box. Reports precision/recall on the test set for several thresholds and writes the model.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import benchmark_gelgenie as bg  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from autogel import __version__  # noqa: E402
from autogel.checker import FEATURES, MODEL_PATH, band_features, fit_logistic, predict  # noqa: E402

SUBSETS = ["quantitation_ladder_gels", "nathan_gels", "external_gels", "stella_gels_for_finetuning", "matthew_gels_2"]


def split_of(sub: str, img: str) -> str:
    folder = os.path.basename(os.path.dirname(img))
    return "test" if sub == "external_gels" or folder.startswith("test") else "train"


def collect(data: str) -> pd.DataFrame:
    rows = []
    for sub in SUBSETS:
        for img, mp in bg.pairs(data, sub):
            try:
                r = bg.evaluate(img, mp)
            except Exception as e:
                print("skip", img, e)
                continue
            ga, m = r["_ga"], r["_mask"]
            gt = [(rg.centroid[0], rg.centroid[1]) for rg in bg.components(m)]
            feats = band_features(ga)
            for k, b in enumerate(ga.bands):
                lane = ga.lanes[b.lane]
                d = lane.dx(b.y)
                hit = any(lane.x0 + d - 3 <= x <= lane.x1 + d + 3 and b.y0 - 3 <= y <= b.y1 + 3 for y, x in gt)
                rows.append({**feats.iloc[k].to_dict(), "real": int(hit), "subset": sub,
                             "image": os.path.basename(img), "split": split_of(sub, img),
                             "true_bands_in_gel": len(gt)})
            print(f"{split_of(sub, img):5s} {sub[:10]:10s} {os.path.basename(img)[:40]:40s} bands {len(ga.bands)}",
                  flush=True)
    return pd.DataFrame(rows)


def report(df: pd.DataFrame, p: np.ndarray, title: str):
    print(f"\n{title}")
    tot_real = df.real.sum()
    print(f"  {'threshold':>9s} {'kept':>6s} {'precision':>9s} {'real kept':>9s} {'false removed':>13s}")
    for t in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
        keep = p >= t
        prec = df.real[keep].mean() if keep.any() else float("nan")
        print(f"  {t:9.1f} {keep.mean():6.1%} {prec:9.1%} {df.real[keep].sum() / tot_real:9.1%} "
              f"{((~keep) & (df.real == 0)).sum() / max(1, (df.real == 0).sum()):13.1%}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("data")
    ap.add_argument("--features", default=os.path.join(HERE, "..", "benchmark_data", "band_checker_features.csv"))
    ap.add_argument("--write-model", action="store_true")
    a = ap.parse_args(argv)
    if os.path.exists(a.features):
        df = pd.read_csv(a.features)
    else:
        df = collect(a.data)
        df.to_csv(a.features, index=False)
    tr, te = df[df.split == "train"], df[df.split == "test"]
    print(f"train {len(tr)} bands ({tr.real.mean():.0%} real), test {len(te)} bands ({te.real.mean():.0%} real)")
    x = tr[FEATURES].to_numpy(float)
    mean, std = np.nanmean(x, 0), np.nanstd(x, 0) + 1e-9
    xs = np.nan_to_num((x - mean) / std)
    coef, b0 = fit_logistic(xs, tr.real.to_numpy(float))
    model = {"features": FEATURES, "mean": mean.tolist(), "std": std.tolist(), "coef": coef.tolist(),
             "intercept": b0}
    for name, part in (("TEST (GelGenie test images + external labs)", te), ("train", tr)):
        report(part, predict(model, part), name)
    for sub, g in te.groupby("subset"):
        report(g, predict(model, g), f"test subset {sub}")
    print("\ncoefficients (standardised):")
    for f, c in sorted(zip(FEATURES, coef), key=lambda t: -abs(t[1])):
        print(f"  {f:18s} {c:+.2f}")
    if a.write_model:
        pt = predict(model, te)
        keep = pt >= 0.2
        model.update({"about": "autogel band checker: P(detected band is real). Logistic regression.",
                      "trained": datetime.date.today().isoformat(), "autogel_version": __version__,
                      "training_data": "GelGenie Zenodo 10.5281/zenodo.14641949 (CC-BY): train + val images",
                      "test_data": "GelGenie test images + external_gels",
                      "test_precision_all": float(te.real.mean()),
                      "test_precision_at_0.2": float(te.real[keep].mean()),
                      "test_real_kept_at_0.2": float(te.real[keep].sum() / te.real.sum())})
        os.makedirs(os.path.dirname(MODEL_PATH), exist_ok=True)
        with open(MODEL_PATH, "w", encoding="utf-8") as f:
            json.dump(model, f, indent=1)
        print("model written:", MODEL_PATH)


if __name__ == "__main__":
    main()

"""Band checker: a small learned model that estimates how likely a detected band is real.

Trained on bands autogel detected in hand-labelled gels (GelGenie dataset, Zenodo 14641949,
CC-BY): a band counts as real when a hand-outlined band lies inside it. The model is a logistic
regression on a few band features - a few kilobytes, milliseconds per gel, no GPU. Training and
evaluation: ``tools/train_band_checker.py``; coefficients: ``models/band_checker.json``.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy import ndimage

MODEL_PATH = os.path.join(os.path.dirname(__file__), "models", "band_checker.json")

FEATURES = ["log_snr_lane", "log_snr_gel", "log_rel_gel", "log_rel_lane", "coverage", "log_width_rel",
            "log_neighbour", "rf", "log_bands_in_lane", "log_peakedness", "saturated"]


def band_features(ga) -> pd.DataFrame:
    """One row per band in ``ga.bands`` (same order)."""
    if not ga.bands:
        return pd.DataFrame(columns=FEATURES)
    heights = np.array([max(ga._peak_height(b), 1e-9) for b in ga.bands])
    vols = np.asarray(ga._volumes(), float)
    sats = np.asarray(ga._saturated(), float)
    gel_noise = max(ga.gel_noise(), 1e-9)
    ref95 = max(float(np.percentile(heights, 95)), 1e-9)
    widths = np.array([b.y1 - b.y0 + 1 for b in ga.bands], float)
    med_w = max(float(np.median(widths)), 1.0)
    lane_noise, lane_max, lane_n = {}, {}, {}
    for i in range(len(ga.lanes)):
        c = ga.lane_profile(i)["corr"]
        lane_noise[i] = 1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) + 1e-9 if len(c) > 10 else gel_noise
    for b, h in zip(ga.bands, heights):
        lane_max[b.lane] = max(lane_max.get(b.lane, 0.0), h)
        lane_n[b.lane] = lane_n.get(b.lane, 0) + 1
    rows = []
    for b, h, v, s, w in zip(ga.bands, heights, vols, sats, widths):
        # signal at the same rows in the neighbouring lanes (streaks across the gel, smears)
        nb = []
        for j in (b.lane - 1, b.lane + 1):
            if 0 <= j < len(ga.lanes) and ga.lanes[j].tier == ga.lanes[b.lane].tier:
                pr = ga.lane_profile(j)
                sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
                if sel.any():
                    nb.append(max(float(ndimage.gaussian_filter1d(pr["corr"], 1.5)[sel].max()), 0.0))
        lane = ga.lanes[b.lane]
        well = lane.y0 + 10
        rf = (b.y - well) / max(1.0, lane.y1 - well)
        rows.append({
            "log_snr_lane": np.log10(h / lane_noise[b.lane]),
            "log_snr_gel": np.log10(h / gel_noise),
            "log_rel_gel": np.log10(h / ref95),
            "log_rel_lane": np.log10(h / lane_max[b.lane]),
            "coverage": ga.coverage(b),
            "log_width_rel": np.log10(w / med_w),
            "log_neighbour": np.log10((np.mean(nb) if nb else 0.0) / h + 1e-3),
            "rf": float(np.clip(rf, -0.2, 1.2)),
            "log_bands_in_lane": np.log10(lane_n[b.lane]),
            "log_peakedness": np.log10(h * w / max(v, 1e-9) + 1e-6),
            "saturated": float(s > 0),
        })
    return pd.DataFrame(rows, columns=FEATURES)


def load_model(path: str = MODEL_PATH) -> dict | None:
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def predict(model: dict, feats: pd.DataFrame) -> np.ndarray:
    """Probability that each band is real."""
    if not len(feats):
        return np.zeros(0)
    x = feats[model["features"]].to_numpy(float)
    x = np.nan_to_num((x - np.asarray(model["mean"])) / np.asarray(model["std"]), nan=0.0, posinf=5, neginf=-5)
    z = x @ np.asarray(model["coef"]) + model["intercept"]
    return 1 / (1 + np.exp(-np.clip(z, -30, 30)))


def fit_logistic(x: np.ndarray, y: np.ndarray, l2: float = 1.0, iters: int = 50) -> tuple[np.ndarray, float]:
    """Logistic regression by Newton's method with an L2 penalty (no extra dependency)."""
    n, d = x.shape
    xb = np.hstack([x, np.ones((n, 1))])
    w = np.zeros(d + 1)
    reg = np.full(d + 1, l2)
    reg[-1] = 0.0                                   # intercept not penalised
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(xb @ w, -30, 30)))
        g = xb.T @ (p - y) + reg * w
        hess = (xb * (p * (1 - p))[:, None]).T @ xb + np.diag(reg)
        step = np.linalg.solve(hess, g)
        w -= step
        if np.abs(step).max() < 1e-7:
            break
    return w[:-1], float(w[-1])

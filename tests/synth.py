"""Synthetic gel images with known ground truth for tests."""
from __future__ import annotations

import numpy as np
import tifffile


def make_gel(path, n_lanes=10, H=1200, W=1000, spacing=80, lane_w=50, x_first=140,
             well_y=60, lanes=None, bg=1500.0, gradient=800.0, noise=60.0, dust=(),
             clip=None, invert=False, rot90=0, seed=0, dtype=np.uint16):
    """lanes: {lane_index: [(y, amplitude, sigma), ...]} (amplitude = peak counts above bg)."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:H, 0:W].astype(float)
    img = bg + gradient * (yy / H) + 200 * (xx / W)
    # wells: small bright rectangles
    for k in range(n_lanes):
        xc = x_first + k * spacing
        img[well_y - 6:well_y + 6, int(xc - lane_w / 2):int(xc + lane_w / 2)] += 1500
    lanes = lanes or {}
    for k, bands in lanes.items():
        xc = x_first + k * spacing
        prof_x = 1 / (1 + np.exp(-(xx - (xc - lane_w / 2)) / 1.5)) / (1 + np.exp((xx - (xc + lane_w / 2)) / 1.5))
        for (y, amp, sig) in bands:
            img += amp * np.exp(-0.5 * ((yy - y) / sig) ** 2) * prof_x
    for (y, x, amp) in dust:
        img += amp * np.exp(-0.5 * (((yy - y) / 2.0) ** 2 + ((xx - x) / 2.0) ** 2))
    img += rng.normal(0, noise, img.shape)
    top = clip if clip is not None else np.iinfo(dtype).max
    img = np.clip(img, 0, top)
    if invert:
        img = top - img
    img = np.rot90(img, rot90)
    tifffile.imwrite(path, img.astype(dtype))
    return path


def ladder(y_top=200, y_bottom=1100, sizes=(10000, 6000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 100),
           amp=8000, sig=4):
    """Bands whose position is exactly linear in log10(size). Returns (bands, size->y)."""
    ls = np.log10(sizes)
    ys = y_top + (ls[0] - ls) / (ls[0] - ls[-1]) * (y_bottom - y_top)
    return [(float(y), amp, sig) for y in ys], dict(zip(sizes, ys))


def y_for_size(size, y_top=200, y_bottom=1100, smax=10000, smin=100):
    return y_top + (np.log10(smax) - np.log10(size)) / (np.log10(smax) - np.log10(smin)) * (y_bottom - y_top)

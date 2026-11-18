"""Image loading, quality checks and orientation."""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage


@dataclass
class GelImage:
    data: np.ndarray            # 2-D float32, bands BRIGHT, lanes vertical, wells on top
    raw_max: np.ndarray         # 2-D bool, pixel at/above saturation level (same geometry)
    path: str
    dtype: str
    saturation_level: float
    warnings: list[str] = field(default_factory=list)
    orientation: dict = field(default_factory=dict)
    # several channels (colour image, multi-page TIFF, added images such as a marker photo):
    # each oriented like ``data``, bands bright. Empty for a plain grey image.
    channels: list[np.ndarray] = field(default_factory=list)
    channel_masks: list[np.ndarray] = field(default_factory=list)
    channel_names: list[str] = field(default_factory=list)
    raw_shape: tuple = ()


def _read(path: str) -> np.ndarray:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tif", ".tiff"):
        import tifffile
        arr = tifffile.imread(path)
    else:
        from skimage import io
        arr = io.imread(path)
    return np.asarray(arr)


def load_channels(path: str) -> tuple[list[np.ndarray], list[str], str]:
    """All channels of an image: grey -> one, colour -> red/green/blue (if they differ),
    multi-page TIFF -> one per page (e.g. two fluorescence channels or exposures)."""
    arr = np.squeeze(_read(path))
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        ch = arr[..., :3]
        if np.allclose(ch[..., 0], ch[..., 1]) and np.allclose(ch[..., 1], ch[..., 2]):
            return [arr[..., 0]], ["grey"], str(arr.dtype)
        return [ch[..., k] for k in range(3)], ["red", "green", "blue"], str(arr.dtype)
    if arr.ndim == 3:
        n = arr.shape[0]
        return [arr[k] for k in range(n)], [f"page {k + 1}" for k in range(n)], str(arr.dtype)
    return [arr], ["grey"], str(arr.dtype)


def _combine(chans: list[np.ndarray], names: list[str]) -> tuple[np.ndarray, list[str]]:
    """The default single image: colour -> channel mean, pages -> page 1 (unchanged behaviour)."""
    if len(chans) == 1:
        return chans[0], []
    if names[0] == "red":
        m = np.mean([c.astype(float) for c in chans], axis=0).astype(chans[0].dtype)
        return m, ["Colour image: analysed as the mean of red/green/blue. A single colour channel "
                   "can be chosen under 'channel'."]
    return chans[0], [f"Multi-page image ({len(chans)} pages): page 1 is analysed. Other pages can be "
                      "chosen under 'channel'."]


def load_raw(path: str) -> tuple[np.ndarray, list[str], str]:
    chans, names, dtype = load_channels(path)
    arr, warnings = _combine(chans, names)
    ext = os.path.splitext(path)[1].lower()
    if ext in (".jpg", ".jpeg"):
        warnings.append("JPEG input: lossy compression and usually auto-contrast - "
                        "use the original TIFF for quantification.")
    if arr.dtype == np.uint8:
        warnings.append("8-bit image: limited dynamic range; prefer the 16-bit TIFF export.")
    return arr, warnings, str(arr.dtype)


def auto_invert(arr: np.ndarray) -> bool:
    """Dark bands on a light background? (bands are the minority signal)"""
    small = _block_mean(arr.astype(float), max(1, min(arr.shape) // 400))
    med, lo, hi = np.median(small), np.percentile(small, 1), np.percentile(small, 99.9)
    return bool((med - lo) > (hi - med))


def orient_channel(raw: np.ndarray, o: dict, invert: bool) -> tuple[np.ndarray, np.ndarray]:
    """One raw channel -> (oriented float32 image with bright bands, saturation mask)."""
    sat, _ = saturation_level(raw)
    a = raw.astype(np.float32)
    if invert:
        a = np.float32(raw.max()) - a
        mask = raw <= (0 if np.issubdtype(raw.dtype, np.integer) else raw.min())
    else:
        mask = raw >= sat
    return apply_orientation(a, o), apply_orientation(mask, o)


def saturation_level(arr: np.ndarray) -> tuple[float, str]:
    """Estimate the clipping level.

    Many imagers clip below the dtype maximum (e.g. 62846 instead of 65535
    after dark-frame subtraction). If a noticeable number of pixels sit at
    the observed maximum, that maximum is taken as the saturation level.
    """
    if np.issubdtype(arr.dtype, np.integer):
        dmax = float(np.iinfo(arr.dtype).max)
    else:
        dmax = float(arr.max())
    obs = float(arr.max())
    n_at_max = int((arr >= obs).sum())
    if obs >= dmax:
        return dmax, f"dtype maximum ({dmax:.0f})"
    if n_at_max >= 5:
        return obs, f"observed maximum {obs:.0f} reached by {n_at_max} px (treated as clip level)"
    return dmax, f"dtype maximum ({dmax:.0f}); image max {obs:.0f}"


# ----------------------------------------------------------------- orientation
def _block_mean(a: np.ndarray, f: int) -> np.ndarray:
    h, w = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    return a[:h, :w].reshape(h // f, f, w // f, f).mean((1, 3))


def _periodicity(profile: np.ndarray) -> float:
    p = profile - ndimage.gaussian_filter1d(profile, len(profile) / 15)
    p = p - p.mean()
    ac = np.correlate(p, p, "full")[len(p) - 1:]
    if ac[0] <= 0:
        return 0.0
    ac = ac / ac[0]
    lo, hi = max(2, len(p) // 60), max(3, len(p) // 5)
    return float(ac[lo:hi].max())


def row_periodicity(img: np.ndarray) -> np.ndarray:
    """Per row: how periodic the horizontal edges are along x (autocorrelation peak, 0..1).

    A row of wells (one slot per lane) scores ~0.9; band rows, gel edges and noise score lower.
    Rows with almost no edge energy get 0."""
    sm = ndimage.gaussian_filter(img.astype(np.float32), 1.5)
    gy = np.abs(ndimage.sobel(sm, 0))
    gy = np.minimum(gy, np.percentile(gy, 99.5))
    h, w = gy.shape
    lo_lag, hi_lag = max(3, w // 60), max(4, w // 4)
    score = np.zeros(h)
    energy = np.zeros(h)
    from scipy import signal as sig
    for y in range(3, h - 3):
        strip = gy[y - 2:y + 3].mean(0)
        q = strip - ndimage.uniform_filter1d(strip, max(5, w // 10))
        q = q - q.mean()
        ac = np.correlate(q, q, "full")[w - 1:]
        energy[y] = float(strip.mean())
        if ac[0] <= 0:
            continue
        ac = ac / ac[0]
        pk, _ = sig.find_peaks(ac[lo_lag:hi_lag])
        if len(pk):
            score[y] = float(ac[lo_lag:hi_lag][pk].max())
    if (energy > 0).any():
        score[energy < 0.25 * np.median(energy[energy > 0])] = 0
    return ndimage.gaussian_filter1d(score, 1)


def lanes_vertical_vote(small: np.ndarray) -> float | None:
    """Fraction of band-like blobs that are wider than tall (lanes vertical) - None if too few.

    Bands are short blobs elongated across the lane. Long straight lines (tray and gel edges,
    wells as a whole) are ignored, they fooled the earlier gradient-energy rule on real photos."""
    from skimage import measure
    th = ndimage.white_tophat(small, size=21)
    thr = th.mean() + 3 * th.std()
    lab = measure.label(th > thr)
    h, w = small.shape
    horiz = vert = 0.0
    for r in measure.regionprops(lab):
        if r.area < 8 or r.area > 0.01 * h * w:
            continue
        y0, x0, y1, x1 = r.bbox
        bw, bh = x1 - x0, y1 - y0
        if max(bw, bh) > 0.25 * min(h, w):          # long line: frame / tray / gel edge
            continue
        ratio = bw / max(bh, 1)
        if ratio >= 1.6:
            horiz += r.area
        elif ratio <= 1 / 1.6:
            vert += r.area
    if horiz + vert < 50:
        return None
    return horiz / (horiz + vert)


UNSURE_SCORE = 0.1   # |vertical_score| below this: warn that the 90° decision is uncertain


def _edge_max(rp: np.ndarray) -> float:
    """Strongest periodic row in the top or bottom 35 % (where a well row would be)."""
    n = len(rp)
    return float(max(rp[: int(n * 0.35)].max(initial=0), rp[int(n * 0.65):].max(initial=0)))


def lane_comb_periodicity(small: np.ndarray) -> tuple[float, float]:
    """Periodicity of the mean profile along x and along y (autocorrelation peak, 0..1).

    Lanes repeat at a regular spacing across the gel, so the profile across the lanes is
    periodic - smeared lanes make this stronger, not weaker."""
    from scipy import signal as sig

    def per(p: np.ndarray) -> float:
        p = p - ndimage.uniform_filter1d(p, max(5, len(p) // 6))
        p = p - p.mean()
        ac = np.correlate(p, p, "full")[len(p) - 1:]
        if ac[0] <= 0:
            return 0.0
        ac = ac / ac[0]
        lo, hi = max(3, len(p) // 60), max(4, len(p) // 4)
        pk, _ = sig.find_peaks(ac[lo:hi])
        return float(ac[lo:hi][pk].max()) if len(pk) else 0.0

    h, w = small.shape
    c = small[h // 10: h - h // 10, w // 10: w - w // 10]
    return per(c.mean(0)), per(c.mean(1))


def detect_orientation(img: np.ndarray, invert: bool | None = None) -> dict:
    """Return dict(invert, rot90_k, flip_lr, angle) to bring lanes vertical, wells on top.

    Lanes vertical vs. horizontal: a weighted sum of band shape, well row, lane comb and edge
    energy (``vertical_cues``). Upside down: conservative, only when a well row is clearly at
    the bottom and not at the top, because most gel photos are saved upright."""
    a = img.astype(float)
    small = _block_mean(a, max(1, min(a.shape) // 400))
    if invert is None:
        # bands are the minority signal: if the median is closer to the max, the background is bright
        med, lo, hi = np.median(small), np.percentile(small, 1), np.percentile(small, 99.9)
        invert = (med - lo) > (hi - med)
    if invert:
        small = small.max() - small
    # Lane direction: weighted evidence from four cues, each positive for "lanes vertical".
    # Band shape alone was fooled by heavily smeared lanes (their streaks look like tall
    # blobs); a periodic well row and the lane comb are not, because both run *across* lanes.
    vote = lanes_vertical_vote(small)
    th = ndimage.white_tophat(small, size=15)
    th = np.minimum(th, np.percentile(th, 99.9))
    ey = float((ndimage.sobel(th, 0) ** 2).sum())
    ex = float((ndimage.sobel(th, 1) ** 2).sum())
    rp0, rp1 = row_periodicity(small), row_periodicity(np.rot90(small, 1))
    comb_x, comb_y = lane_comb_periodicity(small)
    cues = {"band_shape": 0.0 if vote is None else vote - 0.5,
            "well_row": 0.5 * (_edge_max(rp0) - _edge_max(rp1)),
            "lane_comb": comb_x - comb_y,
            "edge_energy": 0.05 * (1.0 if ey >= ex else -1.0)}
    score = float(sum(cues.values()))
    k = 0 if score > 0 else 1
    # upside down only if a row of wells is clearly at the bottom and not at the top
    rp = rp0 if k == 0 else rp1
    n = len(rp)
    top, bottom = rp[: int(n * 0.35)].max(initial=0), rp[int(n * 0.65):].max(initial=0)
    flipped = bottom > 0.6 and bottom > top + 0.15
    if flipped:
        k += 2
    s = np.rot90(small, k)
    angle = estimate_skew(s)
    return {"invert": bool(invert), "rot90_k": int(k % 4), "flip_lr": False,
            "angle": float(angle), "band_shape_vote": None if vote is None else round(float(vote), 3),
            "well_score_top": round(float(top), 3),
            "well_score_bottom": round(float(bottom), 3),
            "vertical_score": round(score, 3),
            "vertical_cues": {c: round(float(v), 3) for c, v in cues.items()}}


SKEW_MIN_GAIN = 0.15   # a tilt must sharpen the lane profile by >= 15 % over 0 degrees


def estimate_skew(small: np.ndarray, max_deg: float = 3.0, step: float = 0.25) -> float:
    """Angle (deg) that makes lanes most vertical (sharpest lane profile).

    0 unless the best angle is clearly better than no rotation: on a real gel a +1 degree
    'tilt' that sharpened the profile by only 3 % moved the lane comb onto the wrong spacing
    (a real 2.5 degree tilt sharpens it by > 150 %)."""
    h = small.shape[0]
    crop = small[h // 10: h - h // 10]
    scores = {}
    for ang in np.arange(-max_deg, max_deg + 1e-9, step):
        r = ndimage.rotate(crop, ang, reshape=False, order=1, mode="nearest")
        prof = r[r.shape[0] // 8: -r.shape[0] // 8 or None].mean(0)
        scores[round(float(ang), 4)] = float(np.var(np.diff(prof)))
    best_a = max(scores, key=scores.get)
    base = scores.get(0.0, 0.0)
    if abs(best_a) < 0.3 or scores[best_a] < (1 + SKEW_MIN_GAIN) * base:
        return 0.0
    return float(best_a)


def apply_orientation(a: np.ndarray, o: dict, order: int = 1) -> np.ndarray:
    out = np.rot90(a, o.get("rot90_k", 0))
    if o.get("flip_lr"):
        out = out[:, ::-1]
    if o.get("flip_ud"):
        out = out[::-1, :]
    ang = o.get("angle", 0.0)
    if ang:
        if out.dtype == bool:
            out = ndimage.rotate(out.astype(np.uint8), ang, reshape=False, order=0, mode="constant") > 0
        else:
            out = ndimage.rotate(out, ang, reshape=False, order=order, mode="nearest")
    return np.ascontiguousarray(out)


def load_gel(path: str, orientation: dict | None = None, invert: bool | None = None) -> GelImage:
    """Load an image and orient it. ``orientation`` (from a saved session) skips auto-detection."""
    arr, warnings, dtype = load_raw(path)
    chans, names, _ = load_channels(path)
    sat, sat_note = saturation_level(arr)
    sat_mask = arr >= sat
    o = dict(orientation) if orientation else detect_orientation(arr, invert)
    if o.get("angle") is None:          # rotated/flipped by the user: re-estimate the skew
        small = _block_mean(arr.astype(float), max(1, min(arr.shape) // 400))
        if o.get("invert"):
            small = small.max() - small
        o["angle"] = float(estimate_skew(apply_orientation(small, {**o, "angle": 0.0})))
    if not orientation and abs(o.get("vertical_score", 1.0)) < UNSURE_SCORE:
        warnings.append("Lane direction uncertain: check that the lanes run top to bottom "
                        "(otherwise click Rotate 90°).")
    a = arr.astype(np.float32)   # exact for <= 24-bit integers, half the RAM of float64
    if o["invert"]:
        a = np.float32(arr.max()) - a
        # for dark-band stains the "saturated" end is 0 (black)
        sat_mask = arr <= (0 if np.issubdtype(arr.dtype, np.integer) else arr.min())
    data = apply_orientation(a, o)
    mask = apply_orientation(sat_mask, o)
    n_sat = int(sat_mask.sum())
    if n_sat:
        warnings.append(f"{n_sat} saturated pixels (level: {sat_note}). "
                        "Bands containing them are flagged and their intensity is underestimated.")
    g = GelImage(data=data, raw_max=mask, path=path, dtype=dtype,
                 saturation_level=sat, warnings=warnings, orientation=o, raw_shape=tuple(arr.shape))
    if len(chans) > 1:
        for c, n in zip(chans, names):
            d, m = orient_channel(c, o, bool(o["invert"]))
            g.channels.append(d)
            g.channel_masks.append(m)
            g.channel_names.append(n)
    return g

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


def _read(path: str) -> np.ndarray:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tif", ".tiff"):
        import tifffile
        arr = tifffile.imread(path)
    else:
        from skimage import io
        arr = io.imread(path)
    return np.asarray(arr)


def load_raw(path: str) -> tuple[np.ndarray, list[str], str]:
    arr = _read(path)
    warnings = []
    arr = np.squeeze(arr)
    if arr.ndim == 3:
        if arr.shape[-1] in (3, 4):
            ch = arr[..., :3].astype(float)
            if np.allclose(ch[..., 0], ch[..., 1]) and np.allclose(ch[..., 1], ch[..., 2]):
                arr = arr[..., 0]
            else:
                warnings.append("RGB image converted to grey (channel mean) - "
                                "colour images are not reliable for quantification.")
                arr = ch.mean(-1).astype(arr.dtype)
        else:
            warnings.append(f"Multi-page image ({arr.shape[0]} pages): using page 1.")
            arr = arr[0]
    ext = os.path.splitext(path)[1].lower()
    if ext in (".jpg", ".jpeg"):
        warnings.append("JPEG input: lossy compression and usually auto-contrast - "
                        "use the original TIFF for quantification.")
    if arr.dtype == np.uint8:
        warnings.append("8-bit image: limited dynamic range; prefer the 16-bit TIFF export.")
    return arr, warnings, str(arr.dtype)


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


def detect_orientation(img: np.ndarray, invert: bool | None = None) -> dict:
    """Return dict(invert, rot90_k, flip_lr, angle) to bring lanes vertical, wells on top."""
    a = img.astype(float)
    small = _block_mean(a, max(1, min(a.shape) // 400))
    if invert is None:
        # bands are the minority signal: if the median is closer to the max, the background is bright
        med, lo, hi = np.median(small), np.percentile(small, 1), np.percentile(small, 99.9)
        invert = (med - lo) > (hi - med)
    if invert:
        small = small.max() - small
    # Bands, wells and tray edges all run ACROSS the lanes, so intensity changes
    # mostly along the run direction: compare gradient energy along both axes
    from scipy import ndimage as ndi
    th = ndi.white_tophat(small, size=15)
    th = np.minimum(th, np.percentile(th, 99.9))
    gy = ndi.sobel(th, 0)
    gx = ndi.sobel(th, 1)
    ey, ex = float((gy ** 2).sum()), float((gx ** 2).sum())
    px, py = ey, ex   # kept for the log
    k = 0 if ey >= ex else 1
    s = np.rot90(small, k)
    # wells on top: signal centre of mass along the run direction lies nearer the wells
    from scipy import ndimage as ndi
    th = ndi.white_tophat(s, size=15)
    rows = np.minimum(th, np.percentile(th, 99.9)).sum(1)
    com = (rows * np.arange(len(rows))).sum() / max(rows.sum(), 1e-9)
    if com > len(rows) / 2:
        k += 2
    s = np.rot90(small, k)
    angle = estimate_skew(s)
    return {"invert": bool(invert), "rot90_k": int(k % 4), "flip_lr": False,
            "angle": float(angle), "grad_energy_run": px, "grad_energy_cross": py}


def estimate_skew(small: np.ndarray, max_deg: float = 3.0, step: float = 0.25) -> float:
    """Angle (deg) that makes lanes most vertical (sharpest lane profile)."""
    best, best_a = -np.inf, 0.0
    h = small.shape[0]
    crop = small[h // 10: h - h // 10]
    for ang in np.arange(-max_deg, max_deg + 1e-9, step):
        r = ndimage.rotate(crop, ang, reshape=False, order=1, mode="nearest")
        prof = r[r.shape[0] // 8: -r.shape[0] // 8 or None].mean(0)
        score = np.var(np.diff(prof))
        if score > best:
            best, best_a = score, float(ang)
    return best_a if abs(best_a) >= 0.3 else 0.0


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
    sat, sat_note = saturation_level(arr)
    sat_mask = arr >= sat
    o = dict(orientation) if orientation else detect_orientation(arr, invert)
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
    return GelImage(data=data, raw_max=mask, path=path, dtype=dtype,
                    saturation_level=sat, warnings=warnings, orientation=o)

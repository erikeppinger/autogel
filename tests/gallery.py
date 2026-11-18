"""Simulated gel gallery with ground truth: different stains, qualities and artefacts.

Each scenario renders a gel image (TIFF, or JPEG for the JPEG case) and knows where every
band is, how big it is and how much material it holds. ``evaluate`` runs the full autogel
pipeline and scores it against that truth.
"""
from __future__ import annotations

import dataclasses
import os

import numpy as np
import tifffile
from scipy import ndimage

DNA_LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 100]
PROT_LADDER = [250, 150, 100, 75, 50, 37, 25, 20, 15, 10]


@dataclasses.dataclass
class Scenario:
    name: str
    description: str
    unit: str = "bp"
    n_lanes: int = 12
    H: int = 1300
    W: int = 1100
    lane_w: int = 52
    spacing: int = 80
    marker_lanes: tuple = (0, 11)          # 0-based
    stain: str = "bright"                  # bright (EtBr/SYBR) | dark (Coomassie/silver)
    rgb: bool = False                      # blue Coomassie-like colour image
    noise: float = 60.0
    bg: float = 1500.0
    gradient: float = 600.0                # top->bottom background ramp
    vignette: float = 0.0                  # fraction of bg lost at the corners
    amp: float = 7000.0                    # peak height of a 50-ng (or 1-unit) band
    band_sigma: float = 5.0
    tilt_deg: float = 0.0
    smile_px: float = 0.0                  # edge lanes lag behind the centre by this much
    dust: int = 0
    scratches: int = 0
    smear_lanes: tuple = ()
    empty_lanes: tuple = ()
    clip: float | None = None              # saturation level
    eight_bit: bool = False
    jpeg_quality: int | None = None
    rot90: int = 0
    seed: int = 1
    tiers: int = 1                         # gel rows (two combs = 2); H is per gel row
    columns: int = 1                       # gels side by side; W is per gel


SCENARIOS = [
    Scenario("agarose_clean", "agarose / EtBr, good exposure"),
    Scenario("agarose_low_signal", "faint bands, high camera noise", amp=900, noise=110),
    Scenario("agarose_uneven_light", "strong vignetting and illumination gradient",
             gradient=4000, vignette=0.6),
    Scenario("agarose_tilted", "gel placed 2.5° crooked on the tray", tilt_deg=2.5),
    Scenario("agarose_smiling", "smiling bands (edge lanes run 35 px slower)", smile_px=35),
    Scenario("agarose_smiling_3_markers", "smiling bands, third marker in the middle lane",
             smile_px=35, marker_lanes=(0, 6, 11)),
    Scenario("agarose_two_combs", "two gels stacked in one image (two combs, own ladders)",
             tiers=2, H=900),
    Scenario("agarose_side_by_side", "two gels side by side on the tray (own ladders), divider between",
             columns=2, n_lanes=10, W=900, marker_lanes=(0, 9)),
    Scenario("agarose_dirty", "dust specks and scratches", dust=120, scratches=8),
    Scenario("agarose_smears_empty", "degraded-DNA smears and two empty lanes",
             smear_lanes=(3, 7), empty_lanes=(5, 9)),
    Scenario("agarose_overexposed", "overexposed: bright bands clipped", amp=26000, clip=30000),
    Scenario("agarose_jpeg_8bit", "8-bit JPEG export (quality 70)", eight_bit=True, jpeg_quality=70),
    Scenario("agarose_rotated", "image saved rotated by 90° and upside down", rot90=1),
    Scenario("agarose_many_narrow", "26 narrow lanes (comb with many teeth)", n_lanes=26,
             W=1500, lane_w=34, spacing=54, marker_lanes=(0, 25)),
    Scenario("agarose_few_wide", "6 wide lanes, one marker", n_lanes=6, W=800, lane_w=90,
             spacing=120, marker_lanes=(0,)),
    Scenario("sds_coomassie", "SDS-PAGE, Coomassie (dark bands on light, blue RGB)", unit="kDa",
             stain="dark", rgb=True, n_lanes=10, W=950, lane_w=60, spacing=85,
             marker_lanes=(0,), band_sigma=7, bg=2000, gradient=300, amp=6000, noise=40),
    Scenario("sds_coomassie_grey", "SDS-PAGE, Coomassie scanned in greyscale, 8-bit", unit="kDa",
             stain="dark", eight_bit=True, n_lanes=10, W=950, lane_w=60, spacing=85,
             marker_lanes=(0, 9), band_sigma=7, bg=2000, amp=6000, noise=40),
]


def _y_for(size, sc: Scenario, y_top=200, y_bot=None):
    lad = PROT_LADDER if sc.unit == "kDa" else DNA_LADDER
    y_bot = y_bot or sc.H - 160
    lmax, lmin = np.log10(lad[0]), np.log10(lad[-1])
    return y_top + (lmax - np.log10(size)) / (lmax - lmin) * (y_bot - y_top)


def truth_lanes(sc: Scenario) -> dict[int, list[tuple[float, float]]]:
    """lane -> [(size, amount)], amount in ng (DNA) or arbitrary units (protein)."""
    rng = np.random.default_rng(sc.seed)
    lad = PROT_LADDER if sc.unit == "kDa" else DNA_LADDER
    out = {}
    for k in range(sc.n_lanes):
        if k in sc.marker_lanes:
            out[k] = [(s, 50.0) for s in lad]
        elif k in sc.empty_lanes:
            out[k] = []
        else:
            n = rng.integers(1, 4)
            lo, hi = np.log10(lad[-2]), np.log10(lad[1])
            sizes = sorted(10 ** rng.uniform(lo, hi, n), reverse=True)
            # keep bands apart (>= 12 % in size) so the truth is unambiguous
            sz = []
            for s in sizes:
                if all(abs(np.log10(s / t)) > 0.05 for t in sz):
                    sz.append(float(s))
            out[k] = [(s, float(rng.uniform(20, 120))) for s in sz]
    return out


def _signal(sc: Scenario, truth: dict) -> np.ndarray:
    H, W = sc.H, sc.W
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    signal = np.zeros((H, W), np.float32)
    x_first = (W - (sc.n_lanes - 1) * sc.spacing) / 2
    for k in range(sc.n_lanes):
        xc = x_first + k * sc.spacing
        e1 = np.clip(-(xx - (xc - sc.lane_w / 2)) / 1.5, -60, 60)
        e2 = np.clip((xx - (xc + sc.lane_w / 2)) / 1.5, -60, 60)
        across = 1 / (1 + np.exp(e1)) / (1 + np.exp(e2))
        smile = sc.smile_px * ((xc - W / 2) / (W / 2)) ** 2          # edges lag behind
        bow = 3.0 * ((xx - xc) / (sc.lane_w / 2)) ** 2                # bands slightly curved
        for size, amount in truth[k]:
            y = _y_for(size, sc) - smile
            sig = sc.band_sigma * (1 + 0.6 * (y - 200) / H)            # lower bands broader
            # integrated signal is proportional to the amount: broader band -> lower peak
            signal += (sc.amp * amount / 50 * sc.band_sigma / sig
                       * np.exp(-0.5 * ((yy - y + bow) / sig) ** 2) * across)
        if k in sc.smear_lanes:
            y0, y1 = _y_for(PROT_LADDER[1] if sc.unit == "kDa" else 8000, sc), _y_for(
                PROT_LADDER[-2] if sc.unit == "kDa" else 300, sc)
            ramp = np.clip((yy - y0) / (y1 - y0), 0, 1) * (yy < y1)
            signal += 0.18 * sc.amp * np.sin(np.pi * ramp) * across
        # well
        signal += 900 * ((yy > 70) & (yy < 82)) * across
    return signal


def render(sc: Scenario, folder: str) -> tuple[str, dict]:
    """Render the scenario; with tiers > 1 several gels are stacked (lanes numbered on)."""
    rng = np.random.default_rng(sc.seed)
    truth, parts = {}, []
    gap = 60                                   # tray divider between gels lying side by side
    for t in range(sc.tiers):
        row = []
        for c in range(sc.columns):
            g = t * sc.columns + c
            sub = dataclasses.replace(sc, seed=sc.seed + g)
            tr = truth_lanes(sub)
            if c:
                div = np.zeros((sc.H, gap), np.float32)
                div[:, gap // 2 - 3: gap // 2 + 3] = 4000.0   # bright vertical tray edge
                row.append(div)
            row.append(_signal(sub, tr))
            truth.update({g * sc.n_lanes + k: v for k, v in tr.items()})
        parts.append(np.hstack(row))
    signal = np.vstack(parts)
    H, W = signal.shape
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    bg = sc.bg + sc.gradient * yy / H
    if sc.vignette:
        r2 = ((xx - W / 2) / (W / 2)) ** 2 + ((yy - H / 2) / (H / 2)) ** 2
        bg = bg * (1 - sc.vignette * np.clip(r2 / 2, 0, 1))
    img = bg + signal
    for _ in range(sc.dust):
        y, x, a = rng.uniform(0, H), rng.uniform(0, W), rng.uniform(2000, 20000)
        r = rng.uniform(1.0, 2.5)
        img += a * np.exp(-0.5 * (((yy - y) / r) ** 2 + ((xx - x) / r) ** 2))
    for _ in range(sc.scratches):
        x0, y0 = rng.uniform(0, W), rng.uniform(0, H)
        ang = rng.uniform(0, np.pi)
        d = np.abs((xx - x0) * np.sin(ang) - (yy - y0) * np.cos(ang))
        along = np.abs((xx - x0) * np.cos(ang) + (yy - y0) * np.sin(ang))
        img += 3000 * np.exp(-0.5 * (d / 0.8) ** 2) * (along < rng.uniform(30, 120))
    img += rng.normal(0, sc.noise, img.shape).astype(np.float32)
    if sc.tilt_deg:
        img = ndimage.rotate(img, sc.tilt_deg, reshape=False, order=1, mode="nearest")
    top = sc.clip if sc.clip else 65535.0
    img = np.clip(img, 0, top)
    if sc.stain == "dark":
        img = top - img
    if sc.eight_bit or sc.jpeg_quality:
        img = np.clip(img / top * 255, 0, 255).astype(np.uint8)
    else:
        img = img.astype(np.uint16)
    if sc.rgb:  # Coomassie blue: weak in blue channel, strong absorption in red/green
        f = img.astype(np.float32) / (255 if img.dtype == np.uint8 else 65535)
        rgb = np.dstack([f ** 1.6, f ** 1.3, f ** 0.6])
        img = (rgb * (255 if img.dtype == np.uint8 else 65535)).astype(img.dtype)
    img = np.rot90(img, sc.rot90 * 3 if sc.rot90 else 0)
    if sc.rot90:
        img = img[::-1]          # also upside down
    os.makedirs(folder, exist_ok=True)
    if sc.jpeg_quality:
        from PIL import Image
        path = os.path.join(folder, sc.name + ".jpg")
        Image.fromarray(np.ascontiguousarray(img)).save(path, quality=sc.jpeg_quality)
    else:
        path = os.path.join(folder, sc.name + ".tif")
        tifffile.imwrite(path, np.ascontiguousarray(img))
    return path, truth


def _orientation_effect(sc: Scenario, o: dict) -> tuple[bool, bool]:
    """Apply render's transform and then the tool's to a coordinate grid: are lanes mirrored
    (left/right) or is the gel upside down relative to the upright rendered gel?"""
    from autogel.imageio import apply_orientation
    Wt = sc.W * sc.columns + 60 * (sc.columns - 1)
    yy, xx = np.mgrid[0:sc.H * sc.tiers, 0:Wt]
    grids = []
    for g in (yy, xx):
        t = np.rot90(g, sc.rot90 * 3 if sc.rot90 else 0)
        if sc.rot90:
            t = t[::-1]
        grids.append(apply_orientation(t.astype(np.float32), {**o, "angle": 0.0}, order=0))
    gy, gx = grids
    if gy.shape != (sc.H * sc.tiers, Wt):
        return False, True
    mirrored = gx[0, 0] > gx[0, -1]
    upside = gy[0, 0] > gy[-1, 0]
    return bool(mirrored), bool(upside)


def evaluate(sc: Scenario, folder: str, export: bool = False) -> dict:
    from autogel.analysis import GelAnalysis, Params
    path, truth = render(sc, folder)
    ga = GelAnalysis(path, Params(n_lanes=sc.n_lanes, n_tiers=sc.tiers, n_columns=sc.columns))
    ga.detect_lanes()
    ga.detect_bands()
    lad = PROT_LADDER if sc.unit == "kDa" else DNA_LADDER
    n_gels = sc.tiers * sc.columns
    markers = [g * sc.n_lanes + m for g in range(n_gels) for m in sc.marker_lanes]
    n_total = sc.n_lanes * n_gels
    ga.set_marker([m + 1 for m in markers], sizes=lad, masses=[50.0] * len(lad), unit=sc.unit)
    df = ga.quantify()
    col = f"size_{sc.unit}"
    res = {"scenario": sc.name, "description": sc.description, "lanes_expected": n_total,
           "lanes_found": len(ga.lanes),
           "warnings": "; ".join(ga.gel.warnings + [w for c in ga.calibrations for w in c.warnings])}
    # ladder bands per marker
    res["ladder_bands"] = "/".join(str(int((df.lane == m + 1).sum())) for m in markers) \
        + f" of {len(lad)}"
    # sample bands: match true bands (lane order as rendered) to detected bands by size
    # where did the rendered (upright) lanes end up after the tool's orientation?
    flip, upside_down = _orientation_effect(sc, ga.gel.orientation)
    res["wells_on_top"] = not upside_down
    size_err, found, missed, extra, ratio_err = [], 0, 0, 0, []
    for k in range(n_total):
        if k in markers or len(ga.lanes) != n_total:
            continue
        kk = (k // sc.n_lanes) * sc.n_lanes + (sc.n_lanes - 1 - k % sc.n_lanes) if flip else k
        sub = df[df.lane == kk + 1]
        det = list(zip(sub[col], sub.volume)) if len(sub) else []
        matched = []
        for size, amount in truth[k]:
            if not det:
                missed += 1
                continue
            j = int(np.argmin([abs(np.log10(d[0] / size)) if np.isfinite(d[0]) else 9 for d in det]))
            e = abs(det[j][0] / size - 1)
            if e < 0.15:
                found += 1
                size_err.append(e)
                matched.append((amount, det[j][1]))
                det.pop(j)
            else:
                missed += 1
        extra += len(det)
        if len(matched) >= 2:      # relative quantification inside the lane
            (a1, v1), (a2, v2) = matched[0], matched[1]
            ratio_err.append(abs((v1 / v2) / (a1 / a2) - 1))
    res.update({"bands_found": found, "bands_missed": missed, "extra_bands": extra,
                "median_size_error_pct": round(100 * float(np.median(size_err)), 1) if size_err else None,
                "max_size_error_pct": round(100 * float(np.max(size_err)), 1) if size_err else None,
                "median_ratio_error_pct": round(100 * float(np.median(ratio_err)), 1) if ratio_err else None,
                "orientation": ga.gel.orientation, "path": path})
    if export:
        from autogel.export import export_all
        res["out"] = export_all(ga, os.path.join(folder, sc.name + "_autogel"))
    return res

"""Core gel analysis: lanes, profiles, background, bands, calibration, quantification.

Coordinates: the oriented image has lanes running top->bottom (wells on top).
``x`` = column, ``y`` = row, both in oriented-image pixels.
"""
from __future__ import annotations

import contextlib
import datetime as _dt
import json
import os
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd
from scipy import ndimage, signal

from . import __version__
from .imageio import GelImage, _combine, auto_invert, load_channels, load_gel, orient_channel, reduce
from .markers import PRESETS, MarkerCalibration, _align, assign_marker

# shoulder splitting (Params.split_shoulders): the smaller part must carry at least this share of the
# larger one (1:5) - real bands are slightly lopsided and a small Gaussian fits their tail; and the
# second curvature minimum must stand out of the curvature noise by this factor
AUTO_SOURCES = ("auto", "ladder")   # bands placed by the program (the peak fit may split them)
LADDER_SPLIT_RSS = 0.5         # ladder-guided split: two peaks must leave < 50 % of one peak's residual
LADDER_SPLIT_MIN_SHARE = 0.03  # ... and the hidden band must hold >= 3 % of the pair's fitted area
LADDER_SPLIT_SHARE_TOL = 2.0   # ... and, with masses, within 2x of the datasheet share
DECONV_MIN_R2 = 0.97            # touching bands: peak fit used only when it explains this share
SHOULDER_MIN_SHARE = 0.2
# curved lanes: a traced bend is applied only above this fraction of the lane width
BEND_MIN_FRAC = 0.15
SHOULDER_D2_SNR = 6.0


def find_gel_region(a: np.ndarray) -> tuple[int, int, int, int] | None:
    """Rectangle of the gel / membrane: the largest region around the image centre whose
    background looks like the centre's. Tray, frame, film edge and table around it differ.
    Returns None when that region is (almost) the whole image."""
    H, W = a.shape
    f = max(1, min(H, W) // 300)
    small = a[: (H // f) * f, : (W // f) * f].reshape(H // f, f, W // f, f).mean((1, 3))
    small = ndimage.median_filter(small, size=5)             # bands are thin: removed by the filter
    h, w = small.shape
    centre = small[h // 4: 3 * h // 4, w // 4: 3 * w // 4]
    ref = float(np.median(centre))
    spread = float(np.median(np.abs(centre - ref))) * 1.4826 + 1e-9
    rng = float(np.percentile(small, 99) - np.percentile(small, 1)) + 1e-9
    tol = max(5 * spread, 0.12 * rng)
    mask = np.abs(small - ref) < tol
    mask = ndimage.binary_closing(mask, iterations=3)
    mask = ndimage.binary_fill_holes(mask)
    lab, n = ndimage.label(mask)
    if n == 0:
        return None
    c = lab[h // 2, w // 2]
    if c == 0:
        sizes = ndimage.sum(mask, lab, range(1, n + 1))
        c = int(np.argmax(sizes)) + 1
    comp = lab == c
    # robust rectangle: rows / columns where most of the line belongs to the region
    rows = np.where(comp.mean(1) > 0.5 * comp.mean(1).max())[0]
    cols = np.where(comp.mean(0) > 0.5 * comp.mean(0).max())[0]
    if not len(rows) or not len(cols):
        return None
    y0, y1, x0, x1 = rows[0] * f, (rows[-1] + 1) * f, cols[0] * f, (cols[-1] + 1) * f
    if (y1 - y0) > 0.95 * H and (x1 - x0) > 0.95 * W:
        return None                                          # gel fills the image: nothing to crop
    if (y1 - y0) < 0.3 * H or (x1 - x0) < 0.3 * W:
        return None                                          # implausibly small: don't trust it
    return int(y0), int(y1), int(x0), int(x1)


# ------------------------------------------------------------------ data model
@dataclass
class Params:
    invert: bool | None = None          # None = auto (bright bands on dark = False)
    n_lanes: int | None = None          # None = auto (per gel row)
    lane_method: str = "auto"           # auto (better of the two) | bands (from the bands) | comb (regular)
    n_tiers: int = 1                    # gel rows = combs stacked in one image
    tier_rows: list[float] | None = None  # split rows between gel rows (None = automatic)
    n_columns: int = 1                  # gels side by side in one image
    column_splits: list[float] | None = None  # split columns between them (None = automatic)
    lane_width_frac: float = 0.8        # integration width as fraction of lane width
    bg_method: str = "rolling_ball"     # rolling_ball | lane_edges | none
    bg_radius: float | None = None      # px; None = image height / 40
    smooth_sigma: float = 1.5           # px, profile smoothing for band detection only
    min_snr: float = 6.0                # band prominence in noise SDs
    min_rel_height: float = 0.01        # band prominence relative to brightest band in lane
    min_gel_rel: float = 0.0            # optional: band peak relative to the gel's strong bands
                                        # (95th pct.). ~0.03 removes noise "bands" in empty lanes of
                                        # noisy/JPEG images, but would also drop faint ladder bands
                                        # on gels with very bright samples - therefore off by default
    # automatic noise floor: drop a band only if it is BOTH weak relative to the gel's strong bands
    # (< noise_floor_rel x 95th pct.) AND close to the gel's noise (< noise_floor_snr x gel noise).
    # Faint ladder bands next to bright samples stay (far above noise); noise "bands" in empty
    # lanes of JPEG/noisy images go. Benchmark: -28 % false bands for -1 % real bands.
    noise_floor: bool = True
    noise_floor_rel: float = 0.07
    noise_floor_snr: float = 10.0
    # band checker (learned, checker.py): every band gets band_probability; optionally bands below
    # checker_filter are removed after detection - never in marker lanes (faint ladder bands
    # can score low next to bright samples)
    checker_filter: float = 0.0
    min_band_distance: int = 6          # px
    skip_well_px: int = 25              # rows below the well ignored by band detection
    bound_frac: float = 0.1             # band edges where signal falls to this fraction of peak
    min_coverage: float = 0.4           # band must cover this fraction of the lane width (rejects dust)
    marker_preset: str | None = None
    marker_sizes: list[float] | None = None
    marker_masses: list[float] | None = None
    marker_unit: str = "bp"
    marker_lanes: list[int] = field(default_factory=list)  # 1-based lane numbers
    # valid size range (resolution range of the gel). None = automatic (whole ladder).
    size_min: float | None = None
    size_max: float | None = None
    # ladder load: datasheet masses refer to marker_ref_load; you loaded marker_load (same unit, e.g. µl)
    marker_load: float | None = None
    marker_ref_load: float | None = None
    sample_loads_ul: dict[str, float] = field(default_factory=dict)  # lane number (str) -> µl loaded
    # mass curve: "size" = response per ng depends on fragment size (measured on the ladder:
    # small and very large fragments give less signal per ng); "linear" = one slope through origin
    mass_model: str = "size"
    # several channels / images (colour image, multi-page TIFF, added image such as a marker photo)
    channel: int | None = None          # analysed channel; None = default (colour mean / page 1),
                                        # -1 = all channels combined (two-colour blots: detect both)
    marker_channel: int | None = None   # marker lanes are read from this channel (e.g. a white-light
                                        # photo of a blot showing the prestained ladder)
    extra_images: list[dict] = field(default_factory=list)  # added images: {"path", "name", "invert"}
    # comparing lanes (see compare.py)
    group_tol: float = 0.04             # same band in two lanes: sizes within 4 % (or Rf within 0.02)
    norm_method: str = "none"           # none | loading_control | total_lane (Western blots)
    lc_size: float | None = None        # loading-control band size, e.g. 42 kDa (actin)
    lc_tol: float = 0.1                 # relative size tolerance for finding it
    lc_channel: int | None = None       # read the loading control / total signal from this channel
    target_channel: int | None = None   # measure the bands to normalise in this channel (two-colour)
    reference_lane: int | None = None   # 1-based lane for fold change (e.g. the untreated control)
    rna_integrity: bool = False         # RNA gels: 28S/18S ratio and degradation per lane
    # overlapping bands (touching integration ranges): split the shared signal by fitted Gaussian
    # peak areas instead of at the valley (a faint band next to a strong one is otherwise skewed)
    deconvolve: bool = True
    # curved lanes (fan-out): estimate per gel how lanes drift sideways with depth and follow it.
    # Opt-in: on the lab's straight Test gels, where lower strips hold only the two ladder lanes,
    # the estimate still found a false drift - so it is off unless the user switches it on.
    curved_lanes: bool = False
    # analysis region: only the gel / membrane, not the tray, frame or film edge around it.
    # crop = [y0, y1, x0, x1] set by the user; auto_crop = find it automatically
    crop: list[float] | None = None
    auto_crop: bool = False
    # unresolved doublets (a shoulder, no valley): split a band in two when the profile shows two
    # curvature minima AND two Gaussians fit much better than one
    split_shoulders: bool = False
    # marker lanes: where the ladder expects a band inside a neighbour's band (e.g. NEB 900 bp in
    # the bright 1000 bp band), fit two peaks there and split (new part: source "ladder")
    ladder_guided: bool = True


@dataclass
class Lane:
    x: float
    width: float
    y0: float
    y1: float
    name: str = ""
    source: str = "auto"
    tier: int = 0                       # gel row (0 = top comb)
    column: int = 0                     # gel side by side (0 = left)
    bend: list[list[float]] | None = None  # curved lane: [[y, dx], ...] centre offset vs. y (None = straight)

    @property
    def x0(self):
        return self.x - self.width / 2

    @property
    def x1(self):
        return self.x + self.width / 2

    def dx(self, y):
        """Sideways offset of the lane centre at row(s) y (0 for a straight lane)."""
        if not self.bend:
            return np.zeros_like(np.asarray(y, float)) if np.ndim(y) else 0.0
        ys, ds = zip(*self.bend)
        out = np.interp(y, ys, ds)
        return out if np.ndim(y) else float(out)

    def x_at(self, y):
        return self.x + self.dx(y)

    def contains(self, x: float, y: float, pad: float = 0.0) -> bool:
        c = self.x_at(y)
        return c - self.width / 2 - pad <= x <= c + self.width / 2 + pad

    def outline(self, n: int = 12) -> np.ndarray:
        """Polygon (row, col) of the lane: left edge top->bottom, right edge bottom->top."""
        ys = np.linspace(self.y0, self.y1, n if self.bend else 2)
        c = self.x_at(ys)
        left = np.c_[ys, c - self.width / 2]
        right = np.c_[ys[::-1], (c + self.width / 2)[::-1]]
        return np.vstack([left, right])


@dataclass
class Band:
    lane: int        # 0-based index into lanes
    y: float         # peak row
    y0: float        # integration start row
    y1: float        # integration end row (inclusive)
    source: str = "auto"


# ---------------------------------------------------------------- lanes from the bands (D31)
# User idea: find clear bands (and wells) across the gel first and use them as lane markers -
# each lane gets the position and width of its own bands, so lanes of different width (narrow
# ladder wells next to wide sample wells) and irregular gaps are handled. The regular comb fit
# stays as the second candidate; the one that explains the band blobs better wins (_lane_fit_score).
BLOB_K = 4.0             # blob threshold: median + 4 robust SDs of the band-feature image inside the gel
BLOB_MIN_SUPPORT = 0.1   # a lane needs >= 10 % of a typical lane's band signal (edge artefacts, specks)


def _extent(prof: np.ndarray) -> tuple[int, int]:
    """[e0, e1) of the gel in an intensity profile; dark margins are trimmed only when they are
    clearly darker than the gel and >= 3 % of the length on that side."""
    n = len(prof)
    sm = ndimage.gaussian_filter1d(np.asarray(prof, float), max(2.0, n / 200))
    lo, hi = np.percentile(sm, 2), np.percentile(sm, 60)
    if hi - lo <= 0:
        return 0, n
    idx = np.where(sm > lo + 0.35 * (hi - lo))[0]
    if len(idx) < 0.2 * n:
        return 0, n
    e0, e1, pad = int(idx[0]), int(idx[-1]) + 1, int(0.01 * n)
    e0 = e0 + pad if e0 >= 0.03 * n else 0
    e1 = e1 - pad if n - e1 >= 0.03 * n else n
    return e0, e1


def _band_blobs(bandimg: np.ndarray, data: np.ndarray) -> list[tuple[float, float, float, float]]:
    """Clear band-like blobs (bands, wells): (x0, x1, y, weight), wider than tall."""
    from skimage import measure
    b = ndimage.gaussian_filter(bandimg, 1.5)               # pixel noise would give speckle
    H, W = b.shape
    r0, r1 = _extent(data.mean(1))
    c0, c1 = _extent(data.mean(0))
    inner = b[r0:r1, c0:c1]
    if inner.size == 0:
        return []
    med = float(np.median(inner))
    thr = med + BLOB_K * 1.4826 * float(np.median(np.abs(inner - med))) + 1e-9
    lab = measure.label(ndimage.binary_opening(b > thr, np.ones((2, 3))))
    out = []
    for r in measure.regionprops(lab, intensity_image=b):
        y0, x0, y1, x1 = r.bbox
        w, h = x1 - x0, y1 - y0
        if w < max(8, W / 150) or w < 1.5 * h or r.area < 0.4 * w * max(2, h / 2):
            continue
        # width at a quarter of the blob's height across its rows (the threshold edge includes the
        # band's soft flanks: a 28 px ladder lane measured 41 px)
        xa, xb = max(0, x0 - 5), min(W, x1 + 5)
        row = bandimg[y0:y1, xa:xb].mean(0)
        base = float(min(row[0], row[-1]))
        on = np.where(row > base + 0.25 * (row.max() - base))[0]
        if len(on):
            x0, x1 = xa + int(on[0]), xa + int(on[-1]) + 1
        out.append((float(x0), float(x1), (y0 + y1) / 2, float(r.intensity_mean) * r.area))
    return out


def _lanes_from_blobs(bl: list, n: int | None = None) -> list[tuple[float, float]]:
    """Lanes (centre, width) from blobs that line up: a blob joins the lane it overlaps by at
    least half its (or the lane's) width. Empty lanes are filled where a gap is a multiple of the
    typical spacing; with ``n`` the n best supported lanes are kept."""
    if not bl:
        return []
    mw = float(np.median([x1 - x0 for x0, x1, _, _ in bl]))
    lanes: list[dict] = []
    for x0, x1, _, w in sorted((b for b in bl if b[1] - b[0] <= 2.5 * mw), key=lambda b: b[0] + b[1]):
        best = None
        for L in lanes:
            lx0, lx1 = np.median(L["x0"]), np.median(L["x1"])
            ov = min(x1, lx1) - max(x0, lx0)
            if ov >= 0.5 * min(x1 - x0, lx1 - lx0) and (best is None or ov > best[0]):
                best = (ov, L)
        if best:
            best[1]["x0"].append(x0)
            best[1]["x1"].append(x1)
            best[1]["s"] += w
            best[1]["n"] += 1
        else:
            lanes.append({"x0": [x0], "x1": [x1], "s": w, "n": 1})
    lanes = [L for L in lanes if L["n"] >= 2] or lanes
    ms = float(np.median([L["s"] for L in lanes]))
    lanes = [L for L in lanes if L["s"] >= BLOB_MIN_SUPPORT * ms]
    res = sorted((float(np.median(L["x0"]) + np.median(L["x1"])) / 2,
                   float(np.median(L["x1"]) - np.median(L["x0"])), L["s"]) for L in lanes)
    if len(res) >= 3:
        sp = float(np.median(np.diff([r[0] for r in res])))
        w0 = float(np.median([r[1] for r in res]))
        filled = [res[0]]
        for a, b in zip(res[:-1], res[1:]):
            k = int(round((b[0] - a[0]) / sp))
            if k >= 2 and abs((b[0] - a[0]) / k - sp) < 0.25 * sp:
                filled += [(a[0] + j * (b[0] - a[0]) / k, w0, 0.0) for j in range(1, k)]
            filled.append(b)
        res = filled
    if n and len(res) > n:
        res = sorted(sorted(res, key=lambda r: -r[2])[:n])
    return [(x, w) for x, w, _ in res]


def _lane_fit_score(lanes: list[tuple[float, float]], bl: list) -> float:
    """How well lanes (centre, width) explain the band blobs: each blob counts with the share of
    its width covered by the lane holding its centre; blobs side by side in the same row of one
    lane (that lane holds two real lanes) count against it, twice."""
    tot = sum(w for *_, w in bl)
    if not lanes or tot <= 0:
        return 0.0
    good = clash = 0.0
    members: dict[int, list] = {}
    for x0, x1, y, w in bl:
        c = (x0 + x1) / 2
        k = next((k for k, (x, wd) in enumerate(lanes) if abs(c - x) <= wd / 2), None)
        if k is None:
            continue
        x, wd = lanes[k]
        good += w * max(0.0, min(x1, x + wd / 2) - max(x0, x - wd / 2)) / max(1.0, x1 - x0)
        members.setdefault(k, []).append((x0, x1, y, w))
    for m in members.values():
        m.sort(key=lambda t: t[2])
        for i, a in enumerate(m):
            for b_ in m[i + 1:]:
                if b_[2] - a[2] > 4:
                    break
                if b_[0] > a[1] or a[0] > b_[1]:
                    clash += min(a[3], b_[3])
    return (good - 2 * clash) / tot


class GelAnalysis:
    """Holds the full state of one gel analysis plus a log of every step."""

    def __init__(self, path: str, params: Params | None = None, orientation: dict | None = None):
        self.params = params or Params()
        self.log: list[dict] = []
        self.listeners: list = []
        self.gel: GelImage = load_gel(path, orientation=orientation, invert=self.params.invert)
        self.path = os.path.abspath(path)
        self._default = (self.gel.data, self.gel.raw_max)      # image analysed when no channel is chosen
        self.lanes: list[Lane] = []
        self.bands: list[Band] = []
        self.well_y: float = 0.0
        self.calibrations: list[MarkerCalibration] = []
        self.mass_fit: dict | None = None
        self._ladder_refined: set[int] = set()     # marker lanes already refined (once per detection)
        self._profiles: dict[int, dict] = {}
        H = self.gel.data.shape[0]
        if self.params.bg_radius is None:
            self.params.bg_radius = max(15.0, H / 40)
        self._log("load", {"path": self.path, "orientation": self.gel.orientation,
                           "dtype": self.gel.dtype, "shape": list(self.gel.data.shape),
                           "saturation_level": self.gel.saturation_level,
                           "channels": self.channel_names()},
                  "; ".join(self.gel.warnings) or "ok")
        for e in self.params.extra_images:                      # replaying a session
            self._attach_image(e)
        if self.params.channel is not None:
            self.set_channel(self.params.channel)

    # ------------------------------------------------------------------ log
    def _log(self, step: str, params: dict | None = None, result: str = ""):
        e = {"time": _dt.datetime.now().isoformat(timespec="seconds"),
             "step": step, "params": params or {}, "result": result}
        self.log.append(e)
        for fn in list(getattr(self, "listeners", [])):   # e.g. the GUI status line
            with contextlib.suppress(Exception):  # a broken listener must never break the analysis
                fn(e)

    @property
    def data(self) -> np.ndarray:
        return self.gel.data

    # ------------------------------------------------------------- channels
    def channel_names(self) -> list[str]:
        return list(self.gel.channel_names) or ["grey"]

    def _channel_arrays(self, ch: int | None):
        if ch is None or not self.gel.channels:
            return self._default
        if ch == -1:                                 # all channels: each scaled to its strong signal
            if getattr(self, "_composite", None) is None:
                parts = []
                for c in self.gel.channels:
                    lo, hi = np.percentile(c[::4, ::4], [5, 99.9])
                    parts.append((c - lo) / max(float(hi - lo), 1e-9))
                comp = np.max(parts, axis=0).astype(np.float32) * 10000
                mask = np.logical_or.reduce(self.gel.channel_masks)
                self._composite = (comp, mask)
            return self._composite
        return self.gel.channels[ch], self.gel.channel_masks[ch]

    def set_channel(self, ch: int | None):
        """Analyse another channel (lanes and bands are kept; profiles are re-read)."""
        if ch is not None and not (0 <= ch < len(self.gel.channels) or (ch == -1 and len(self.gel.channels) > 1)):
            raise ValueError(f"channel {ch} does not exist (image has {len(self.gel.channels)} channels)")
        self.params.channel = ch
        self.gel.data, self.gel.raw_max = self._channel_arrays(ch)
        self._profiles = {}
        self._log("set_channel", {"channel": ch},
                  "default image" if ch is None else "all channels combined" if ch == -1
                  else f"channel {ch + 1}: {self.gel.channel_names[ch]}")

    def add_image(self, path: str, name: str | None = None, invert: bool | None = None) -> int:
        """Add a second image of the same gel/membrane as an extra channel, e.g. the white-light
        photo showing a prestained marker next to a chemiluminescence image, or a total-protein
        stain. It must have the same size; the gel's orientation is applied to it."""
        e = {"path": os.path.abspath(path), "name": name or os.path.splitext(os.path.basename(path))[0],
             "invert": invert}
        k = self._attach_image(e)
        self.params.extra_images.append(e)
        self._log("add_image", dict(e), f"channel {k + 1}: {self.gel.channel_names[k]}")
        return k

    def _attach_image(self, e: dict) -> int:
        chans, names, _ = load_channels(e["path"])
        raw, _ = _combine(chans, names)
        if tuple(raw.shape) != tuple(self.gel.raw_shape):
            raise ValueError(f"{os.path.basename(e['path'])} is {raw.shape[1]}x{raw.shape[0]} px but the gel image "
                             f"is {self.gel.raw_shape[1]}x{self.gel.raw_shape[0]} px - both must come from the "
                             "same camera position without cropping.")
        raw = reduce(raw, self.gel.reduced)            # same reduction as the gel image
        inv = auto_invert(raw) if e.get("invert") is None else bool(e["invert"])
        e["invert"] = inv
        if not self.gel.channels:            # plain grey gel so far: it becomes channel 1
            self.gel.channels.append(self._default[0])
            self.gel.channel_masks.append(self._default[1])
            self.gel.channel_names.append("image")
        d, m = orient_channel(raw, self.gel.orientation, inv)
        self.gel.channels.append(d)
        self.gel.channel_masks.append(m)
        self.gel.channel_names.append(e["name"])
        self._profiles = {}
        self._composite = None
        return len(self.gel.channels) - 1

    def set_marker_channel(self, ch: int | None):
        """Read the marker lanes from another channel; their bands are detected there again."""
        if ch is not None and not 0 <= ch < len(self.gel.channels):
            raise ValueError(f"channel {ch} does not exist")
        self.params.marker_channel = ch
        self._profiles = {}
        self._log("set_marker_channel", {"marker_channel": ch},
                  "same image as the samples" if ch is None else f"channel {ch + 1}: {self.gel.channel_names[ch]}")
        mk = [L - 1 for L in self.params.marker_lanes if 0 < L <= len(self.lanes)]
        if mk and self.lanes:
            self.detect_bands(lanes=mk)

    def _lane_arrays(self, i: int, channel="auto"):
        if channel != "auto":
            return self._channel_arrays(channel)
        mc = self.params.marker_channel
        if mc is not None and (i + 1) in self.params.marker_lanes and mc < len(self.gel.channels):
            return self._channel_arrays(mc)
        return self.gel.data, self.gel.raw_max

    # ------------------------------------------------------------------ lanes
    def detect_lanes(self) -> list[Lane]:
        """Find lanes. Several gels in one image - side by side (columns) and/or stacked
        (gel rows, i.e. two combs) - are handled separately: each gel gets its own wells, lane
        comb and lane names (A1.., B1.. in reading order: left to right, top to bottom)."""
        a = self.data
        H, W = a.shape
        # Band-like features only: a vertical top-hat keeps short horizontal
        # structures (bands, wells) and removes long vertical ones (tray edges,
        # smears, illumination gradients).
        bandimg = self._band_feature_image()
        cy0, cy1, cx0, cx1 = self.region()
        columns = [(max(ca, cx0), min(cb, cx1)) for ca, cb in self.column_bounds(bandimg)]
        tiers = [(max(ra, cy0), min(rb, cy1)) for ra, rb in self.tier_bounds()]
        columns = [c for c in columns if c[1] - c[0] > 20]
        tiers = [t for t in tiers if t[1] - t[0] > 20]
        n_gels = len(columns) * len(tiers)
        lanes: list[Lane] = []
        info = []
        for t, (ra, rb) in enumerate(tiers):
            for col, (ca, cb) in enumerate(columns):
                h, w = rb - ra, cb - ca
                prof = bandimg[ra + int(h * .02): rb - int(h * .02), ca:cb].mean(0)
                # baseline window must be wider than any lane (wide lanes were split in two otherwise)
                prof = prof - ndimage.minimum_filter1d(ndimage.gaussian_filter1d(prof, 2), max(5, w // 4))
                comb = prof / (prof.max() + 1e-12)
                peaks, spacing = self._fit_lane_grid(comb, self.params.n_lanes)
                width = 0.8 * spacing
                cand = [(float(p), width) for p in peaks]
                how = "regular comb"
                if self.params.lane_method in ("auto", "bands"):
                    bl = _band_blobs(bandimg[ra:rb, ca:cb], a[ra:rb, ca:cb])
                    from_bands = _lanes_from_blobs(bl, self.params.n_lanes)
                    if self.params.n_lanes and len(from_bands) != int(self.params.n_lanes):
                        from_bands = []                    # the user's lane count is binding
                    if from_bands and (self.params.lane_method == "bands" or
                                       _lane_fit_score(from_bands, bl) > _lane_fit_score(cand, bl) + 0.02):
                        cand, how = from_bands, "from the bands"
                peaks = [x for x, _ in cand]
                well = ra + self._find_well(a[ra:rb, ca:cb], peaks, float(np.median([w for _, w in cand]))
                                            if cand else width)
                y0 = max(float(ra), well - 10)
                y1 = float(rb - 1 - (int(H * .005) if t == len(tiers) - 1 and rb == H else 2))
                prefix = chr(65 + t * len(columns) + col) if n_gels > 1 else ""
                lanes += [Lane(x=float(ca + x), width=float(w), y0=y0, y1=y1, name=f"{prefix}{i + 1}",
                               tier=t, column=col) for i, (x, w) in enumerate(cand)]
                ws = [w for _, w in cand]
                sp_ = float(np.median(np.diff(peaks))) if len(peaks) > 1 else spacing
                info.append(f"{prefix + ': ' if prefix else ''}{len(cand)} lanes {how}, width "
                            + (f"{min(ws):.0f}-{max(ws):.0f}" if ws and max(ws) - min(ws) >= 2 else f"{width:.0f}")
                            + f" px, spacing {sp_:.1f} px, well row {well:.0f}")
        self.lanes = lanes
        n_bent = self.trace_lanes(bandimg) if self.params.curved_lanes else 0
        if n_bent:
            info.append(f"{n_bent} curved lane(s) traced")
        self.well_y = lanes[0].y0 + 10 if lanes else 0.0
        self.bands = []
        self._profiles = {}
        self._log("detect_lanes", {"n_lanes_param": self.params.n_lanes, "n_tiers": len(tiers),
                                   "tier_bounds": [[int(x), int(y)] for x, y in tiers],
                                   "n_columns": len(columns),
                                   "column_bounds": [[int(x), int(y)] for x, y in columns]},
                  "; ".join(info))
        return self.lanes

    def trace_lanes(self, bandimg: np.ndarray | None = None, strips: int = 8) -> int:
        """Curved / fanned-out lanes, estimated per gel from ALL its lanes at once.

        In each horizontal strip the column profile of the gel is matched to the profile of the
        whole gel height, allowing a shift and a stretch about the gel centre. Lanes that fan out
        (or a photo taken at an angle) show up as a stretch that changes with depth. Shift and
        stretch are smoothed over the strips (quadratic in y) and turned into a bend per lane.
        Single-lane tracing was tried first: noise, dust and vignetting gave false bends.
        A bend is applied only when it is consistent over the strips and larger than
        BEND_MIN_FRAC of the lane width. Returns the number of curved lanes."""
        b = self._band_feature_image() if bandimg is None else bandimg
        H, W = b.shape
        for l in self.lanes:
            l.bend = None
        n = 0
        groups: dict[tuple, list[int]] = {}
        for i, l in enumerate(self.lanes):
            groups.setdefault((l.tier, l.column), []).append(i)
        for idx in groups.values():
            if len(idx) < 3:
                continue
            lanes = [self.lanes[i] for i in idx]
            xs = np.array([l.x for l in lanes])
            sp = float(np.median(np.diff(np.sort(xs))))
            ca, cb = max(0, int(xs.min() - sp)), min(W, int(xs.max() + sp) + 1)
            top = int(max(l.y0 for l in lanes) + 10 + self.params.skip_well_px)
            bot = int(min(l.y1 for l in lanes))
            if bot - top < 8 * strips or cb - ca < 20:
                continue
            xc = float(xs.mean())
            xx = np.arange(ca, cb, dtype=float)

            def norm(p):
                p = p - p.mean()
                return p / (np.linalg.norm(p) + 1e-12)
            tmpl = b[top:bot, ca:cb].mean(0)
            scales = np.arange(0.95, 1.0501, 0.0025)
            max_shift = 0.35 * sp
            shifts = np.arange(-max_shift, max_shift + 1e-9, 1.0)
            edges = np.linspace(top, bot, strips + 1)
            pts = []
            for k in range(strips):
                r0, r1 = int(edges[k]), int(edges[k + 1])
                prof = norm(b[r0:r1, ca:cb].mean(0))
                best = (-2.0, 0.0, 1.0)
                for c in scales:
                    for s in shifts:
                        # strip position x corresponds to template position xc + (x - xc - s) / c
                        t = np.interp(xc + (xx - xc - s) / c, xx, tmpl)
                        r = float(prof @ norm(t))
                        if r > best[0]:
                            best = (r, s, c)
                if best[0] > 0.5:                      # the strip shows the lane pattern
                    pts.append((0.5 * (r0 + r1), best[1], best[2], best[0]))
            if len(pts) < 5:
                continue
            y, s_, c_, w = map(np.asarray, zip(*pts))
            half = float(xs.max() - xc)
            # offsets of the outermost lane: shift + stretch x half span
            d_out = s_ + (c_ - 1) * half
            cf_s = np.polyfit(y, s_, 2, w=w)
            cf_c = np.polyfit(y, c_, 2, w=w)
            res = d_out - (np.polyval(cf_s, y) + (np.polyval(cf_c, y) - 1) * half)
            ys = np.linspace(top, bot, 9)
            model = np.polyval(cf_s, ys) + (np.polyval(cf_c, ys) - 1) * half
            rms = float(np.sqrt(np.mean(res ** 2)))
            if rms > 0.1 * sp or np.ptp(model) < 4 * max(rms, 0.5):
                continue
            for l in lanes:
                ny = np.linspace(l.y0, l.y1, 9)
                yy = np.clip(ny, top, bot)
                d = np.polyval(cf_s, yy) + (np.polyval(cf_c, yy) - 1) * (l.x - xc)
                if np.ptp(d) > BEND_MIN_FRAC * l.width:
                    l.bend = [[float(a), float(e)] for a, e in zip(ny, d)]
                    n += 1
        if n:
            self._profiles = {}
        return n

    def region(self) -> tuple[int, int, int, int]:
        """(y0, y1, x0, x1) of the part of the image that is analysed: the user's crop, the
        automatically found gel/membrane (Params.auto_crop) or the whole image."""
        H, W = self.data.shape
        p = self.params
        if p.crop:
            y0, y1, x0, x1 = (int(round(v)) for v in p.crop)
        elif p.auto_crop:
            if getattr(self, "_auto_region", None) is None:
                self._auto_region = find_gel_region(self.data)
            y0, y1, x0, x1 = self._auto_region or (0, H, 0, W)
        else:
            return 0, H, 0, W
        y0, y1 = max(0, min(y0, y1)), min(H, max(y0, y1))
        x0, x1 = max(0, min(x0, x1)), min(W, max(x0, x1))
        if y1 - y0 < 40 or x1 - x0 < 40:
            return 0, H, 0, W
        return y0, y1, x0, x1

    def set_crop(self, crop: list[float] | None, reason: str = "user"):
        self.params.crop = [float(v) for v in crop] if crop else None
        msg = ("whole image" if not crop else
               "analysed region rows {:.0f}-{:.0f}, columns {:.0f}-{:.0f}".format(*crop))
        self._log("set_crop", {"crop": self.params.crop, "reason": reason}, msg)

    def column_bounds(self, bandimg: np.ndarray | None = None) -> list[tuple[int, int]]:
        """Column ranges of gels lying side by side. Manual splits win; otherwise searched."""
        W = self.data.shape[1]
        n = max(1, int(self.params.n_columns or 1))
        if n == 1:
            return [(0, W)]
        cols = sorted(self.params.column_splits or [])
        if len(cols) != n - 1:
            cols = self.find_column_splits(n, bandimg)
            self.params.column_splits = cols
        edges = [0] + [int(c) for c in cols] + [W]
        return [(edges[i], edges[i + 1]) for i in range(n) if edges[i + 1] - edges[i] > 20]

    def find_column_splits(self, n: int, bandimg: np.ndarray | None = None) -> list[float]:
        """Split between gels side by side: the widest stripes without band signal in the middle
        part of the image (tray dividers and gel edges are long vertical lines, already removed
        by the band filter, so they do not count as signal)."""
        b = self._band_feature_image() if bandimg is None else bandimg
        W = b.shape[1]
        prof = ndimage.gaussian_filter1d(b.mean(0), max(2.0, W / 200))
        thr = 0.1 * float(np.percentile(prof, 95))
        low = prof < thr
        runs, start = [], None
        for x in range(W + 1):
            if x < W and low[x]:
                start = x if start is None else start
            elif start is not None:
                if start > 0.15 * W and x < 0.85 * W:          # gaps at the image border don't count
                    runs.append((x - start, (start + x - 1) / 2))
                start = None
        runs.sort(reverse=True)
        cuts = sorted(c for _, c in runs[: n - 1])
        if len(cuts) < n - 1:                                  # no clear gap: equal parts
            cuts = [W * (k + 1) / n for k in range(n - 1)]
        return [float(c) for c in cuts]

    def set_column_splits(self, cols: list[float], reason: str = "user"):
        self.params.column_splits = sorted(float(c) for c in cols)
        self.params.n_columns = len(cols) + 1
        self._log("set_columns", {"columns": self.params.column_splits, "reason": reason},
                  f"{self.params.n_columns} gels side by side, split at column(s) "
                  + ", ".join(f"{c:.0f}" for c in self.params.column_splits))

    def _find_well(self, a: np.ndarray, peaks, width) -> float:
        """Well row inside one gel row: every lane has a well, only some have bands -> median
        over lanes, then the TOPMOST clear peak in the top 20 % (not the strongest one)."""
        H = a.shape[0]
        # below the gel's top edge: where the image shows dark space above the gel, the gel's
        # bright rim right after the rise would be taken for the wells
        head = a[: max(10, int(H * .3))]
        cols = [head[:, max(0, int(p - width / 2)): int(p + width / 2) + 1].mean(1) for p in peaks]
        prof = ndimage.gaussian_filter1d(np.median(cols, axis=0) if cols else head.mean(1), 2)
        lo, mid = float(prof[: max(3, int(H * .01))].mean()), float(np.median(prof))
        g0 = 0
        if lo < 0.3 * mid:                           # dark space above the gel (lab gel: 6 %)
            rise = int(np.argmax(prof > lo + 0.5 * (mid - lo)))
            g0 = rise + int(0.015 * H)
        top = a[g0: g0 + max(5, int(H * .2))]
        per_lane = [top[:, max(0, int(p - width / 2)): int(p + width / 2) + 1].mean(1) for p in peaks]
        rowp = np.median(per_lane, axis=0) if per_lane else top.mean(1)
        rowp = ndimage.gaussian_filter1d(rowp, 2)
        pk, pr = signal.find_peaks(rowp, prominence=0)
        noise = 3 * np.std(np.diff(rowp))
        good = []
        if len(pk):
            thr = max(noise, 0.15 * pr["prominences"].max())
            good = [k for k, pp in zip(pk, pr["prominences"]) if pp > thr]
        return float(g0 + good[0]) if good else float(g0 + int(H * .01))

    def tier_bounds(self) -> list[tuple[int, int]]:
        """Row ranges of the gel rows. Manual split rows win; otherwise they are searched."""
        H = self.data.shape[0]
        n = max(1, int(self.params.n_tiers or 1))
        if n == 1:
            return [(0, H)]
        rows = sorted(self.params.tier_rows or [])
        if len(rows) != n - 1:
            rows = self.find_tier_splits(n)
            self.params.tier_rows = rows
        edges = [0] + [int(r) for r in rows] + [H]
        return [(edges[i], edges[i + 1]) for i in range(n) if edges[i + 1] - edges[i] > 20]

    def find_tier_splits(self, n: int) -> list[float]:
        """Rows where a lower comb's wells start.

        Wells are a row of slots across the whole gel, one per lane, so their horizontal edges
        form a strongly *periodic* pattern along x. For every row we measure that periodicity
        and pick the n-1 strongest rows (well apart) in the lower 80 % of the image. The split
        is placed a little above the well so the well belongs to the lower gel."""
        a = self.data
        H, W = a.shape
        f = max(1, min(H, W) // 600)
        sm = ndimage.gaussian_filter(a[::f, ::f].astype(np.float32), 1.5)
        gy = np.abs(ndimage.sobel(sm, 0))
        gy = np.minimum(gy, np.percentile(gy, 99.5))
        h, w = gy.shape
        lo_lag, hi_lag = max(3, w // 60), max(4, w // 4)
        score = np.zeros(h)
        energy = np.zeros(h)
        for y in range(3, h - 3):
            strip = gy[y - 2:y + 3].mean(0)
            q = strip - ndimage.uniform_filter1d(strip, max(5, w // 10))
            q = q - q.mean()
            ac = np.correlate(q, q, "full")[w - 1:]
            if ac[0] <= 0:
                continue
            ac = ac / ac[0]
            pk, _ = signal.find_peaks(ac[lo_lag:hi_lag])
            energy[y] = float(strip.mean())
            if len(pk):
                # periodicity only: well edges repeat in every lane (~0.9), band rows don't
                # (bright bands must not win just because they have strong edges)
                score[y] = float(ac[lo_lag:hi_lag][pk].max())
        score[energy < 0.25 * np.median(energy[energy > 0])] = 0
        score = ndimage.gaussian_filter1d(score, 1)
        cand = []
        min_sep = h / (n + 1)
        y_lo = int(h * 0.2)
        order = np.argsort(score[y_lo:])[::-1] + y_lo
        for y in order:
            if all(abs(y - c) > min_sep for c in cand):
                cand.append(int(y))
            if len(cand) == n - 1:
                break
        return sorted(float(c * f - 25) for c in cand)

    def set_tier_rows(self, rows: list[float], reason: str = "user"):
        self.params.tier_rows = sorted(float(r) for r in rows)
        self.params.n_tiers = len(rows) + 1
        self._log("set_tiers", {"rows": self.params.tier_rows, "reason": reason},
                  f"{self.params.n_tiers} gel rows, split at row(s) "
                  + ", ".join(f"{r:.0f}" for r in self.params.tier_rows))

    def _band_feature_image(self) -> np.ndarray:
        a = self.data
        size = int(max(31, a.shape[0] / 40)) | 1
        th = ndimage.white_tophat(a, size=(size, 1))
        # suppress single bright specks (dust) and cap saturated bands
        th = ndimage.median_filter(th, size=(3, 3))
        return np.minimum(th, np.percentile(th, 99.5))

    @staticmethod
    def _fit_comb_n(c: np.ndarray, n: int) -> tuple[list[float], float]:
        """Fit a comb with exactly ``n`` equidistant teeth (the lane count the user confirmed).

        Searches every spacing and offset: score = signal on the teeth minus signal halfway
        between them. Works with empty lanes and smears (they keep their comb position),
        where spacing estimates from the visible lanes alone fail. Each tooth is then moved
        onto its lane's plateau if a clear one is nearby (lanes are not perfectly equidistant).
        """
        W = len(c)
        x = np.arange(W)
        if n == 1:
            cs = ndimage.gaussian_filter1d(c, 3)
            p = int(np.argmax(cs))
            above = np.where(cs > 0.5 * cs[p])[0]
            width = float(max(8, above.max() - above.min())) if len(above) else W / 4
            return [float(p)], width / 0.8
        # every lane must count equally: cap the profile at the typical lane height, otherwise
        # two teeth inside one very bright lane (ladder) can beat the true comb
        c_raw = c
        # local normalisation: divide by a sliding maximum (floored at 20 % of the global max so
        # empty background stays low). Bright and dim halves of a gel then weigh the same, and
        # a tight comb squeezed into the bright half no longer wins.
        env = ndimage.maximum_filter1d(ndimage.gaussian_filter1d(c, 2), max(9, W // 8))
        c = c / np.maximum(env, 0.2 * float(c.max()) + 1e-12)
        # the gel may fill only part of the frame, so small spacings must be allowed ...
        s_lo = max(6.0, W / 150)
        s_hi = (W - 1) / (n - 1)
        # ... but then a narrow comb squeezed into part of the gel can out-score the true one.
        # Coverage term: the comb must span the region where the lanes' signal is.
        # signal region = 10th-90th percentile of the profile's mass (sparse noise, specks and
        # tray edges outside the gel hardly shift it)
        mass = np.cumsum(np.clip(c_raw, 0, None))
        if mass[-1] > 0:
            sig_span = float(np.searchsorted(mass, 0.9 * mass[-1]) - np.searchsorted(mass, 0.1 * mass[-1]))
        else:
            sig_span = float(W)
        teeth = np.arange(n)
        best = (-np.inf, s_lo, 0.0)
        cache: dict[float, np.ndarray] = {}
        for sp in np.arange(s_lo, s_hi, 0.25 if s_hi - s_lo < 150 else 0.5):
            sig = round(max(1.0, sp / 8) * 2) / 2
            if sig not in cache:
                cache[sig] = ndimage.gaussian_filter1d(c, sig)
            cs = cache[sig]
            x0 = np.arange(0, W - (n - 1) * sp, 1.0)
            if not len(x0):
                continue
            pos = x0[:, None] + sp * teeth[None, :]
            mid = pos[:, :-1] + sp / 2
            coverage = min(1.0, (n - 1) * sp / max(sig_span * 0.85, 1.0))
            score = (np.interp(pos, x, cs).mean(1) - np.interp(mid, x, cs).mean(1)) * coverage
            k = int(np.argmax(score))
            if score[k] > best[0]:
                best = (float(score[k]), float(sp), float(x0[k]))
        _, s, x0 = best
        cs = ndimage.gaussian_filter1d(c_raw, max(1.0, s / 8))
        out = []
        for k in range(n):
            p = x0 + k * s
            lo, hi = int(max(0, p - s / 3)), int(min(W, p + s / 3 + 1))
            seg = cs[lo:hi]
            if len(seg) < 3 or seg.max() - seg.min() < 0.02:
                out.append(float(p))            # empty lane: keep the comb position
                continue
            m = lo + int(np.argmax(seg))
            half = seg.min() + 0.5 * (cs[m] - seg.min())
            left = right = m
            while left > lo and cs[left - 1] > half:
                left -= 1
            while right < hi - 1 and cs[right + 1] > half:
                right += 1
            out.append(float(0.5 * (left + right)))
        return out, float(s)

    @staticmethod
    def _fit_lane_grid(c: np.ndarray, n_lanes: int | None = None,
                       max_gap: int = 3) -> tuple[list[float], float]:
        """Fit equidistant lanes (comb geometry) to the lane profile ``c``.

        Score of a grid = summed signal on lane positions minus signal halfway
        between them; this penalises half, double and odd-multiple spacings.
        Returns refined lane centres and the spacing.
        """
        W = len(c)
        x = np.arange(W)
        smin, smax = max(8.0, W / 80), W / 4
        if n_lanes:
            return GelAnalysis._fit_comb_n(c, int(n_lanes))

        def grid_score(sp: float) -> float:
            cs_ = ndimage.gaussian_filter1d(c, sp / 8)
            best_ = -np.inf
            for ph in np.arange(0, sp, max(1.0, sp / 24)):
                pos = np.arange(ph, W, sp)
                mid = pos + sp / 2
                n = int((mid < W).sum())
                if n == 0:
                    continue
                # mean contrast x sqrt(#positions): mean alone favours huge spacings that
                # hit only the brightest lanes, a plain sum favours tiny ones (noise)
                d = np.interp(pos[:n], x, cs_) - np.interp(mid[:n], x, cs_)
                best_ = max(best_, float(d.mean() * np.sqrt(n)))
            return best_

        # Spacing: first autocorrelation peak that is at least half as prominent
        # as the strongest one (later peaks are multiples of the true spacing).
        q = c - c.mean()
        ac = np.correlate(q, q, "full")[W - 1:]
        ac = ac / (ac[0] + 1e-12)
        lo, hi = int(smin), int(max(smin + 2, smax))
        pk_, pr_ = signal.find_peaks(ac[lo:hi], prominence=0.01)
        if len(pk_):
            prom_ = pr_["prominences"]
            first = int(np.where(prom_ >= 0.5 * prom_.max())[0][0])
            s0 = float(pk_[first] + lo)
        else:
            s0 = float(W / 10)
        # local refinement on the comb-grid score
        fine = np.arange(max(smin, s0 * 0.95), s0 * 1.05, 0.25)
        s = float(fine[int(np.argmax([grid_score(sp) for sp in fine]))]) if len(fine) else s0
        # The grid gives the spacing; actual centres come from the lane plateaus,
        # because spacing is not perfectly constant across a whole image.
        cs = ndimage.gaussian_filter1d(c, s / 8)
        pk, pr = signal.find_peaks(cs, distance=max(3, int(0.6 * s)), prominence=0.02)
        if len(pk) == 0:
            return [], float(s)
        prom = pr["prominences"]
        centres = []
        for p in pk:
            lo, hi = int(max(0, p - s / 2)), int(min(W, p + s / 2 + 1))
            base = cs[lo:hi].min()
            half = base + 0.5 * (cs[p] - base)
            left = p
            while left > lo and cs[left - 1] > half:
                left -= 1
            right = p
            while right < hi - 1 and cs[right + 1] > half:
                right += 1
            centres.append(0.5 * (left + right))
        centres = np.asarray(centres)

        def fill(kept):
            out = [float(kept[0])]
            for a_, b_ in zip(kept[:-1], kept[1:]):
                n = int(round((b_ - a_) / s))
                for j in range(1, n):
                    out.append(float(a_ + j * (b_ - a_) / n))
                out.append(float(b_))
            return out

        if n_lanes:
            # best run of peaks that yields exactly n_lanes lanes after gap filling,
            # scored by total prominence (ladders/strong lanes anchor the run)
            best_run, best_score = None, -np.inf
            for i in range(len(centres)):
                for j in range(i, len(centres)):
                    nl = len(fill(centres[i:j + 1]))
                    if nl > n_lanes:
                        break
                    if nl == n_lanes:
                        sc = prom[i:j + 1].sum()
                        if sc > best_score:
                            best_run, best_score = (i, j), sc
            if best_run is not None:
                return fill(centres[best_run[0]:best_run[1] + 1]), float(s)
            # fewer visible lanes than requested (empty edge lanes): extend the run
            # by one spacing at a time, alternating right/left, while it fits the image
            out = fill(centres)
            right = True
            while len(out) < n_lanes:
                nxt_r, nxt_l = out[-1] + s, out[0] - s
                fits_r, fits_l = nxt_r <= W - s / 2, nxt_l >= s / 2
                if not (fits_r or fits_l):
                    break
                if (right and fits_r) or not fits_l:
                    out.append(nxt_r)
                else:
                    out.insert(0, nxt_l)
                right = not right
            return out[:n_lanes], float(s)
        else:
            # relative to the median lane: bright ladders must not hide weak sample lanes
            strong = np.where(prom > 0.3 * np.median(prom))[0]
            clusters, cur = [], []
            for k in strong:
                if cur and (centres[k] - centres[cur[-1]]) > (max_gap + 0.5) * s:
                    clusters.append(cur)
                    cur = []
                cur.append(k)
            if cur:
                clusters.append(cur)
            main = max(clusters, key=lambda cl: prom[cl].sum())
            sel = list(range(main[0], main[-1] + 1))
        return fill(centres[sel]), float(s)

    def set_lanes(self, lanes: list[Lane], reason: str = "user edit"):
        old = len(self.lanes)
        self._prev_lanes = list(self.lanes)
        lanes = sorted(lanes, key=lambda l: l.x)
        names = [l.name for l in self.lanes]
        for i, l in enumerate(lanes):
            if not l.name or l.name.isdigit():
                l.name = names[i] if i < len(names) and names[i] and not names[i].isdigit() else str(i + 1)
        self.lanes = lanes
        self._profiles = {}
        # re-attach bands to the lane that now contains them
        old_lanes = getattr(self, "_prev_lanes", None) or []
        new_bands = []
        for b in self.bands:
            if b.lane < len(old_lanes):
                x = old_lanes[b.lane].x_at(b.y)
                hit = [j for j, l in enumerate(lanes) if l.contains(x, b.y)]
                if hit:
                    b.lane = hit[0]
                    new_bands.append(b)
        self.bands = new_bands
        self._log("set_lanes", {"reason": reason,
                                "lanes": [[round(l.x, 1), round(l.width, 1), round(l.y0), round(l.y1)]
                                          + ([round(l.dx(l.y1) - l.dx(l.y0), 1)] if l.bend else []) for l in lanes]},
                  f"{old} -> {len(lanes)} lanes, {len(new_bands)} bands kept")

    def set_lane_names(self, names: list[str]):
        for l, n in zip(self.lanes, names):
            l.name = n
        self._log("set_lane_names", {"names": names})

    # ------------------------------------------------------------- profiles
    def lane_profile(self, i: int, channel="auto") -> dict:
        """Background-corrected profile of lane i. ``channel``: "auto" = the analysed image (the
        marker channel for marker lanes, if set), or a channel index / None for the default image."""
        key = i if channel == "auto" else (i, channel)
        if key in self._profiles:
            return self._profiles[key]
        lane = self.lanes[i]
        a, satmask = self._lane_arrays(i, channel)
        H, W = a.shape
        hw = lane.width * self.params.lane_width_frac / 2
        c0, c1 = max(0, int(round(lane.x - hw))), min(W, int(round(lane.x + hw)) + 1)
        ya, yb = sorted((lane.y0, lane.y1))
        r0, r1 = max(0, int(ya)), min(H, int(yb) + 1)
        r1 = max(r0, r1)
        rows = np.arange(r0, r1)
        off = np.rint(lane.dx(rows.astype(float))).astype(int) if lane.bend else np.zeros(len(rows), int)
        raw = np.zeros(len(rows))
        sat = np.zeros(len(rows), dtype=np.int64)
        for o in np.unique(off):
            m = off == o
            ca, cb = max(0, c0 + o), min(W, c1 + o)
            if cb > ca:
                raw[m] = a[rows[m], ca:cb].sum(1, dtype=np.float64)
                sat[m] = satmask[rows[m], ca:cb].sum(1)
        base = self._baseline(raw, i, (r0, r1, c0, c1), off, a)
        prof = {"y": rows.astype(float), "raw": raw, "baseline": base,
                "corr": raw - base, "sat": sat, "cols": (c0, c1), "offset": off}
        self._profiles[key] = prof
        return prof

    def smear_profile(self, i: int) -> tuple[np.ndarray, np.ndarray]:
        """(rows, signal) of lane i with the background taken from BESIDE the lane. A rolling
        ball treats a broad smear (degraded RNA/DNA) as background and removes it."""
        pr = self.lane_profile(i)
        y = pr["y"]
        if not len(y):
            return y, y
        a = self._lane_arrays(i)[0]
        r0, r1 = int(y[0]), int(y[-1]) + 1
        c0, c1 = pr["cols"]
        base = self._baseline(pr["raw"], i, (r0, r1, c0, c1), pr["offset"], a, method="lane_edges")
        return y, np.clip(pr["raw"] - base, 0, None)

    def _baseline(self, raw: np.ndarray, i: int, box, off: np.ndarray | None = None,
                  a: np.ndarray | None = None, method: str | None = None) -> np.ndarray:
        m = method or self.params.bg_method
        if m == "none":
            return np.zeros_like(raw)
        if m == "lane_edges":
            r0, r1, c0, c1 = box
            a = self.data if a is None else a
            g = max(2, int(self.lanes[i].width * 0.1))
            o = int(np.median(off)) if off is not None and len(off) else 0
            left = a[r0:r1, max(0, c0 + o - 2 * g): max(0, c0 + o - g)]
            right = a[r0:r1, c1 + o + g: c1 + o + 2 * g]
            parts = [p for p in (left, right) if p.size]
            if not parts:
                return np.zeros_like(raw)
            per_px = np.median(np.hstack(parts), axis=1)
            per_px = ndimage.median_filter(per_px, size=max(3, int(self.params.bg_radius) | 1))
            return np.minimum(per_px * (c1 - c0), raw)
        # rolling ball (1-D): morphological opening + smoothing, never above the signal
        size = max(3, int(2 * self.params.bg_radius) | 1)
        sm = ndimage.gaussian_filter1d(raw, 1.0)
        opened = ndimage.grey_opening(sm, size=size)
        base = ndimage.uniform_filter1d(opened, size=max(3, size // 2))
        return np.minimum(base, raw)

    # ---------------------------------------------------------------- bands
    def detect_bands(self, lanes: list[int] | None = None) -> list[Band]:
        idx = range(len(self.lanes)) if lanes is None else lanes
        self.bands = [b for b in self.bands if b.lane not in set(idx)]
        self._ladder_refined -= set(idx)
        total = rejected = 0
        new: list[tuple[Band, float]] = []
        for i in idx:
            pr = self.lane_profile(i)
            corr = pr["corr"]
            y = pr["y"]
            sm = ndimage.gaussian_filter1d(corr, self.params.smooth_sigma)
            noise = 1.4826 * np.median(np.abs(np.diff(corr))) / np.sqrt(2) + 1e-9
            valid = y >= self.lanes[i].y0 + 10 + self.params.skip_well_px
            if not valid.any():
                continue
            thr = max(self.params.min_snr * noise, self.params.min_rel_height * sm[valid].max())
            pk, props = signal.find_peaks(sm, prominence=thr, distance=self.params.min_band_distance,
                                          width=1.5)
            pk = [p for p in pk if valid[p]]
            for j, p in enumerate(pk):
                lo, hi = self._bounds(sm, p, pk[j - 1] if j > 0 else None,
                                      pk[j + 1] if j + 1 < len(pk) else None)
                b = Band(lane=i, y=float(y[p]), y0=float(y[lo]), y1=float(y[hi]))
                if self.coverage(b) < self.params.min_coverage:
                    rejected += 1
                    continue
                new.append((b, float(sm[p])))
        # gel-wide floor: compare with the strong bands of the whole gel (incl. bands kept from
        # other lanes) - a lane-relative threshold cannot judge a lane that holds only noise
        floor = 0
        if new and self.params.noise_floor:
            ref = float(np.percentile([h for _, h in new] + [self._peak_height(b) for b in self.bands], 95))
            gel_noise = self.gel_noise()
            keep = [(b, h) for b, h in new
                    if not (h < self.params.noise_floor_rel * ref and h < self.params.noise_floor_snr * gel_noise)]
            floor = len(new) - len(keep)
            new = keep
        weak = 0
        if new and self.params.min_gel_rel > 0:
            others = [self._peak_height(b) for b in self.bands]
            ref = float(np.percentile([h for _, h in new] + others, 95))
            keep = [(b, h) for b, h in new if h >= self.params.min_gel_rel * ref]
            weak = len(new) - len(keep)
            new = keep
        unlikely = 0
        if new and self.params.checker_filter > 0:
            from .checker import band_features, load_model, predict
            model = load_model()
            if model is not None:
                markers = {L - 1 for L in self.params.marker_lanes}
                kept_before = list(self.bands)
                self.bands = kept_before + [b for b, _ in new]
                p = predict(model, band_features(self))[len(kept_before):]
                self.bands = kept_before
                keep = [(b, h) for (b, h), pi in zip(new, p)
                        if b.lane in markers or pi >= self.params.checker_filter]
                unlikely = len(new) - len(keep)
                new = keep
        split = 0
        if self.params.split_shoulders:
            out = []
            for b, h in new:
                parts = self._split_shoulder(b)
                split += len(parts) - 1
                out += [(p_, h) for p_ in parts]
            new = out
        for b, _ in new:
            self.bands.append(b)
            total += 1
        self.bands.sort(key=lambda b: (b.lane, b.y))
        self._log("detect_bands", {"lanes": None if lanes is None else [l + 1 for l in lanes],
                                   "bg_method": self.params.bg_method, "bg_radius": self.params.bg_radius,
                                   "min_snr": self.params.min_snr, "min_rel_height": self.params.min_rel_height,
                                   "min_gel_rel": self.params.min_gel_rel,
                                   "noise_floor": self.params.noise_floor,
                                   "split_shoulders": self.params.split_shoulders,
                                   "checker_filter": self.params.checker_filter,
                                   "skip_well_px": self.params.skip_well_px},
                  f"{total} bands ({rejected} rejected as narrower than "
                  f"{self.params.min_coverage:.0%} of lane width: dust/specks"
                  + (f"; {weak} below {self.params.min_gel_rel:.0%} of the gel's strong bands" if weak else "")
                  + (f"; {floor} at the noise floor" if floor else "")
                  + (f"; {unlikely} unlikely (band checker)" if unlikely else "")
                  + (f"; {split} shoulder(s) split into two bands" if split else "")
                  + ")")
        return self.bands

    def gel_noise(self) -> float:
        """Typical noise of a lane profile in this gel (median over lanes of the robust
        point-to-point noise), in profile units."""
        vals = []
        for i in range(len(self.lanes)):
            c = self.lane_profile(i)["corr"]
            if len(c) > 10:
                vals.append(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2))
        return float(np.median(vals)) if vals else 0.0

    def _split_shoulder(self, b: Band) -> list[Band]:
        """[b] or two bands if b is an unresolved doublet (see Params.split_shoulders)."""
        from scipy.optimize import curve_fit
        pr = self.lane_profile(b.lane)
        sel = (pr["y"] >= b.y0 - 2) & (pr["y"] <= b.y1 + 2)
        yy, cc = pr["y"][sel], pr["corr"][sel]
        if len(yy) < 12 or pr["sat"][sel].sum() > 0:     # clipped flat top looks like two peaks
            return [b]
        d2 = -ndimage.gaussian_filter1d(cc, 2.0, order=2)
        noise_d2 = 1.4826 * np.median(np.abs(-ndimage.gaussian_filter1d(pr["corr"], 2.0, order=2))) + 1e-9
        pk, _ = signal.find_peaks(d2, prominence=SHOULDER_D2_SNR * noise_d2, distance=4)
        if len(pk) < 2:
            return [b]
        pk = sorted(pk, key=lambda i: -d2[i])[:2]
        g1 = lambda y, a, m, s_: a * np.exp(-0.5 * ((y - m) / s_) ** 2)  # noqa: E731
        g2 = lambda y, a1, m1, s1, a2, m2, s2: g1(y, a1, m1, s1) + g1(y, a2, m2, s2)  # noqa: E731
        amax, w = float(cc.max()), max(1.0, (b.y1 - b.y0) / 4)
        try:
            with np.errstate(all="ignore"):
                p1, _ = curve_fit(g1, yy, cc, p0=[amax, b.y, w], maxfev=3000)
                ya, yb = sorted(float(yy[i]) for i in pk)
                p2, _ = curve_fit(g2, yy, cc, p0=[amax, ya, w / 2, amax / 2, yb, w / 2],
                                  bounds=([0, yy[0], 0.5, 0, yy[0], 0.5], [3 * amax, yy[-1], len(yy), 3 * amax,
                                                                           yy[-1], len(yy)]), maxfev=4000)
        except (RuntimeError, ValueError):
            return [b]
        rss1 = float(((cc - g1(yy, *p1)) ** 2).sum())
        rss2 = float(((cc - g2(yy, *p2)) ** 2).sum())
        a1, m1, s1, a2, m2, s2 = p2
        if m1 > m2:
            a1, m1, s1, a2, m2, s2 = a2, m2, s2, a1, m1, s1
        areas = (a1 * s1, a2 * s2)
        sep_ok = (m2 - m1) >= 1.2 * (s1 + s2) / 2 and (m2 - m1) >= self.params.min_band_distance
        if not (rss2 < 0.4 * rss1 and sep_ok and min(areas) >= SHOULDER_MIN_SHARE * max(areas)):
            return [b]
        # boundary where the two fitted components are equal
        yy_f = np.linspace(m1, m2, 200)
        cut = float(yy_f[np.argmin(np.abs(g1(yy_f, a1, m1, s1) - g1(yy_f, a2, m2, s2)))])
        return [Band(b.lane, float(m1), b.y0, float(np.floor(cut)), b.source),
                Band(b.lane, float(m2), float(np.floor(cut)) + 1, b.y1, b.source)]

    def _peak_height(self, b: Band) -> float:
        """Smoothed background-corrected profile value at the band's peak row."""
        pr = self.lane_profile(b.lane)
        if not len(pr["y"]):
            return 0.0
        sm = ndimage.gaussian_filter1d(pr["corr"], self.params.smooth_sigma)
        return float(sm[int(np.clip(round(b.y - pr["y"][0]), 0, len(sm) - 1))])

    def coverage(self, b: Band) -> float:
        """Fraction of the lane's columns that carry the band (bands span the lane, dust does not)."""
        pr = self.lane_profile(b.lane)
        c0, c1 = pr["cols"]
        o = int(round(self.lanes[b.lane].dx(b.y)))
        img = self._lane_arrays(b.lane)[0]
        W = img.shape[1]
        c0, c1 = max(0, c0 + o), min(W, c1 + o)
        y0, y1 = int(b.y0), int(b.y1) + 1
        block = img[y0:y1, c0:c1].sum(0, dtype=np.float64)
        sel = (pr["y"] >= y0) & (pr["y"] < y1)
        # per-column share of the lane background (independent of neighbouring bands)
        col = block - pr["baseline"][sel].sum() / max(1, c1 - c0)
        if col.size == 0:
            return 0.0
        col = ndimage.uniform_filter1d(col, 3)
        if col.max() <= 0:
            return 0.0
        return float((col > 0.3 * col.max()).mean())

    def _bounds(self, sm, p, prev_p, next_p) -> tuple[int, int]:
        h = sm[p]
        lim = self.params.bound_frac * h
        lo_stop = prev_p + int(np.argmin(sm[prev_p:p + 1])) if prev_p is not None else 0
        hi_stop = p + int(np.argmin(sm[p:next_p + 1])) if next_p is not None else len(sm) - 1
        lo = p
        while lo > lo_stop and sm[lo - 1] > lim:
            lo -= 1
        hi = p
        while hi < hi_stop and sm[hi + 1] > lim:
            hi += 1
        return lo, hi

    def band_for(self, lane: int, y0: float, y1: float, source="user") -> Band:
        """Build a band from integration limits; peak = max of corrected profile inside."""
        pr = self.lane_profile(lane)
        sel = (pr["y"] >= y0) & (pr["y"] <= y1)
        if sel.any():
            y = float(pr["y"][sel][np.argmax(pr["corr"][sel])])
        else:
            y = (y0 + y1) / 2
        return Band(lane=lane, y=y, y0=float(y0), y1=float(y1), source=source)

    def band_at(self, lane: int, y: float, source="user") -> Band:
        """Create a band around row ``y`` using automatic edge finding."""
        pr = self.lane_profile(lane)
        if len(pr["y"]) == 0:   # lane outside the image or of zero height
            return Band(lane=lane, y=float(y), y0=float(y), y1=float(y), source=source)
        sm = ndimage.gaussian_filter1d(pr["corr"], self.params.smooth_sigma)
        p = int(np.clip(round(y - pr["y"][0]), 0, len(sm) - 1))
        w = 4
        p = max(0, p - w) + int(np.argmax(sm[max(0, p - w): p + w + 1]))
        lo, hi = self._bounds(sm, p, None, None)
        return Band(lane=lane, y=float(pr["y"][p]), y0=float(pr["y"][lo]), y1=float(pr["y"][hi]), source=source)

    def set_bands(self, bands: list[Band], reason: str = "user edit"):
        old = {(b.lane, round(b.y0), round(b.y1)) for b in self.bands}
        new = {(b.lane, round(b.y0), round(b.y1)) for b in bands}
        added, removed = new - old, old - new
        self.bands = sorted(bands, key=lambda b: (b.lane, b.y))
        self._log("set_bands", {"reason": reason,
                                "added": sorted([[l + 1, a, b] for l, a, b in added]),
                                "removed": sorted([[l + 1, a, b] for l, a, b in removed])},
                  f"{len(added)} added/changed, {len(removed)} removed/changed, {len(bands)} total")

    # ---------------------------------------------------------- calibration
    def set_marker(self, lanes_1based: list[int], preset: str | None = None,
                   sizes: list[float] | None = None, masses: list[float] | None = None,
                   unit: str | None = None, load: float | None = None, ref_load: float | None = None):
        p = self.params
        before = (p.marker_preset, p.marker_sizes, p.marker_masses, p.marker_unit,
                  list(p.marker_lanes), p.marker_load, p.marker_ref_load)
        if load is not None:
            p.marker_load = load
        if ref_load is not None:
            p.marker_ref_load = ref_load
        if preset and sizes is None:
            pre = PRESETS[preset]
            p.marker_preset, p.marker_sizes = preset, list(pre["sizes"])
            p.marker_masses = list(pre["masses"]) if pre.get("masses") else None
            p.marker_unit = pre["unit"]
        if sizes is not None:
            p.marker_sizes = list(sizes)
            p.marker_masses = list(masses) if masses else None
            p.marker_preset = preset or "custom"
        if unit:
            p.marker_unit = unit
        p.marker_lanes = list(lanes_1based)
        after = (p.marker_preset, p.marker_sizes, p.marker_masses, p.marker_unit,
                 list(p.marker_lanes), p.marker_load, p.marker_ref_load)
        if after == before:
            return
        self._ladder_refined = set()
        if p.marker_channel is not None and list(p.marker_lanes) != before[4]:
            self._profiles = {}
            mk = [L - 1 for L in p.marker_lanes if 0 < L <= len(self.lanes)]
            if mk:
                self.detect_bands(lanes=mk)
        self._log("set_marker", {"lanes": p.marker_lanes, "preset": p.marker_preset,
                                 "sizes": p.marker_sizes, "masses": p.marker_masses, "unit": p.marker_unit,
                                 "load": p.marker_load, "ref_load": p.marker_ref_load})

    def size_range(self) -> tuple[float, float]:
        p = self.params
        # 2 % tolerance: a "10 kb" limit must include NEB's 10002 bp band
        return (float(p.size_min) * 0.98 if p.size_min else 0.0,
                float(p.size_max) * 1.02 if p.size_max else float("inf"))

    def set_size_range(self, lo: float | None, hi: float | None, reason: str = ""):
        if lo and hi and lo > hi:
            lo, hi = hi, lo
        p = self.params
        if (p.size_min, p.size_max) == (lo or None, hi or None):
            return
        p.size_min, p.size_max = lo or None, hi or None
        self._log("set_size_range", {"min": p.size_min, "max": p.size_max, "reason": reason},
                  "automatic (whole ladder)" if not (lo or hi) else f"{lo or 0:g} - {hi or 'inf'}")

    def refine_marker_bands(self) -> list[str]:
        """Ladder-guided split of merged marker bands (D30), once per marker lane after detection.

        The lane's bands are paired with the ladder (:func:`markers._align`). A ladder size
        without band between two paired neighbours gets its row interpolated in log(size). If that
        row lies inside a neighbour's band (or on its flank, just outside its range), the band may
        hold both (a bright band with a faint
        neighbour on its shoulder, e.g. NEB 1000 + 900 bp). Two Gaussians are fitted there, one
        near the band's peak and one near the expected row (+-4 px). If they leave less than
        ``LADDER_SPLIT_RSS`` of one Gaussian's residual and the hidden part holds at least
        ``LADDER_SPLIT_MIN_SHARE`` of the area (with masses: within ``LADDER_SPLIT_SHARE_TOL`` of
        the datasheet's share), the band is split where the two components are
        equal; the new part gets source "ladder". Volumes then come from the usual peak fit or
        valley split of touching bands. Returns what was done (also logged)."""
        from scipy.optimize import curve_fit
        p = self.params
        done: list[str] = []
        if not (p.ladder_guided and p.marker_sizes and p.marker_lanes) or len(p.marker_sizes) < 5:
            return done
        lg = np.log10(np.asarray(p.marker_sizes, float))

        def g1(y, a, m, s_):
            return a * np.exp(-0.5 * ((y - m) / s_) ** 2)

        def g2(y, a1, m1, s1, a2, m2, s2):
            return g1(y, a1, m1, s1) + g1(y, a2, m2, s2)
        for L in p.marker_lanes:
            i = L - 1
            if not 0 <= i < len(self.lanes) or i in self._ladder_refined:
                continue
            self._ladder_refined.add(i)
            lane_b = sorted([b for b in self.bands if b.lane == i], key=lambda b: b.y)
            if len(lane_b) < 5:
                continue
            vol = dict(zip(map(id, self.bands), self._volumes()))
            pairs = _align([b.y for b in lane_b], p.marker_sizes, p.marker_masses, [vol[id(b)] for b in lane_b])
            if len(pairs) < 5:
                continue
            by_size = {j: lane_b[k] for k, j in pairs}
            pr = self.lane_profile(i)
            js = sorted(by_size)
            for a, c in zip(js, js[1:]):
                if c - a != 2:                       # exactly one size missing between two pairs
                    continue
                j, ba, bc = a + 1, by_size[a], by_size[c]
                y_hat = ba.y + (bc.y - ba.y) * (lg[a] - lg[j]) / (lg[a] - lg[c])
                host = ba if ba.y0 <= y_hat <= ba.y1 else bc if bc.y0 <= y_hat <= bc.y1 else None
                lo, hi = (host.y0, host.y1) if host is not None else (0.0, 0.0)
                if host is None:
                    # on the flank, just outside the nearer neighbour's range: fit band + flank
                    near = ba if abs(y_hat - ba.y) <= abs(bc.y - y_hat) else bc
                    edge = near.y1 if near is ba else near.y0
                    if abs(y_hat - edge) < 0.5 * abs(y_hat - near.y):
                        host, reach = near, abs(y_hat - near.y)
                        other = bc if near is ba else ba
                        if near is ba:
                            lo, hi = host.y0, min(y_hat + 0.5 * reach, other.y0 - 1)
                        else:
                            lo, hi = max(y_hat - 0.5 * reach, other.y1 + 1), host.y1
                if host is None or host.source not in AUTO_SOURCES or abs(host.y - y_hat) < 3:
                    continue
                sel = (pr["y"] >= lo) & (pr["y"] <= hi)
                yy, cc = pr["y"][sel], pr["corr"][sel]
                if len(yy) < 8 or cc.max() <= 0:
                    continue
                amax, w = float(cc.max()), max(1.5, (host.y1 - host.y0) / 6)
                try:
                    with np.errstate(all="ignore"):
                        q1, _ = curve_fit(g1, yy, cc, p0=[amax, host.y, w], maxfev=3000)
                        q2, _ = curve_fit(g2, yy, cc, p0=[amax, host.y, w, amax / 5, y_hat, w],
                                          bounds=([0, host.y - 3, 0.8, 0, y_hat - 4, 0.8],
                                                  [3 * amax, host.y + 3, 4 * w, amax, y_hat + 4, 4 * w]),
                                          maxfev=5000)
                except (RuntimeError, ValueError):
                    continue
                rss1 = float(((cc - g1(yy, *q1)) ** 2).sum())
                rss2 = float(((cc - g2(yy, *q2)) ** 2).sum())
                a1, m1, s1, a2, m2, s2 = q2
                share = a2 * s2 / max(a1 * s1 + a2 * s2, 1e-12)
                if not (rss2 < LADDER_SPLIT_RSS * rss1 and share >= LADDER_SPLIT_MIN_SHARE and abs(m2 - m1) >= 3):
                    continue
                if p.marker_masses:                  # the split must roughly agree with the datasheet
                    hj = a if host is ba else c
                    expect = p.marker_masses[j] / max(p.marker_masses[j] + p.marker_masses[hj], 1e-9)
                    if not LADDER_SPLIT_SHARE_TOL ** -1 <= share / max(expect, 1e-9) <= LADDER_SPLIT_SHARE_TOL:
                        continue
                yf = np.linspace(m1, m2, 200)
                cut = float(np.floor(yf[np.argmin(np.abs(g1(yf, a1, m1, s1) - g1(yf, a2, m2, s2)))]))
                if m2 > m1:
                    parts = [Band(i, host.y, lo, cut, host.source), Band(i, float(m2), cut + 1, float(hi), "ladder")]
                else:
                    parts = [Band(i, float(m2), float(lo), cut, "ladder"), Band(i, host.y, cut + 1, hi, host.source)]
                if min(b.y1 - b.y0 for b in parts) < 2:
                    continue
                self.bands.remove(host)
                self.bands += parts
                done.append(f"lane {L}: {p.marker_sizes[j]:g} found inside the band at row {host.y:.0f} "
                            f"(now row {m2:.0f}, {share:.0%} of the pair)")
        if done:
            self.bands.sort(key=lambda b: (b.lane, b.y))
            self._log("refine_marker_bands", {"marker_lanes": p.marker_lanes, "preset": p.marker_preset},
                      "; ".join(done))
        return done

    def calibrate(self) -> list[MarkerCalibration]:
        p = self.params
        self.calibrations, self.mass_fit = [], None
        if not p.marker_sizes or not p.marker_lanes:
            return []
        self.refine_marker_bands()
        vols = self._volumes()
        for L in p.marker_lanes:
            i = L - 1
            if not 0 <= i < len(self.lanes):
                continue
            bl = [(b, v) for b, v in zip(self.bands, vols) if b.lane == i]
            ys = [b.y for b, _ in bl]
            ys_a, sizes, masses, w = assign_marker(ys, p.marker_sizes, p.marker_masses,
                                                   volumes=[v for _, v in bl])
            vmap = {b.y: v for b, v in bl}
            # bands are assigned with the whole ladder; only those inside the valid size
            # range are used for the size curve and the mass curve
            lo, hi = self.size_range()
            keep = [k for k, s in enumerate(sizes) if lo <= s <= hi]
            if len(keep) < 2 and len(sizes) >= 2:
                w = w + [f"Marker lane {L}: fewer than 2 ladder bands inside the size range "
                         f"{lo:g}-{hi:g}; range ignored for this lane."]
                keep = list(range(len(sizes)))
            cal = MarkerCalibration(lane=i, ys=[ys_a[k] for k in keep], sizes=[sizes[k] for k in keep],
                                    masses=[masses[k] for k in keep] if masses else None,
                                    volumes=[vmap[ys_a[k]] for k in keep], warnings=w)
            self.calibrations.append(cal)
        self._fit_mass()
        warn = [w for c in self.calibrations for w in c.warnings]
        self._log("calibrate", {"marker_lanes": p.marker_lanes},
                  ("; ".join(warn) if warn else "ok") +
                  (f"; mass curve ({self.mass_fit['model']}): typical error "
                   f"{self.mass_fit['typical_error']:.0%} (leave-one-out), R²={self.mass_fit['r2']:.3f}, "
                   f"n={self.mass_fit['n']}" if self.mass_fit else ""))
        return self.calibrations

    def mass_scale(self) -> float:
        """Loaded ladder / datasheet reference load. The reference defaults to the
        library entry's reference_load (e.g. NEB 1 kb Plus: 1 µg)."""
        p = self.params
        ref = p.marker_ref_load
        if not ref:
            ref = (PRESETS.get(p.marker_preset or "", {}).get("reference_load") or {}).get("amount")
        if p.marker_load and ref:
            return float(p.marker_load) / float(ref)
        return 1.0

    def _fit_mass(self):
        """Mass standard curve from the unsaturated ladder bands.

        ``size`` model (default): log10(volume / ng) as a smooth function (polynomial in log size,
        degree 1-3 by the number of distinct sizes) - on real gels the signal per ng peaks at
        1.5-3 kb and drops to ~20 % at 100 bp, so one slope mis-estimates small fragments 3-5x.
        ``linear``: volume = slope x ng through the origin.
        The leave-one-out error (typical factor between estimated and datasheet mass) is reported."""
        xs, ys, ss = [], [], []
        k = self.mass_scale()
        sat = self._saturated()
        for c in self.calibrations:
            if not c.masses:
                continue
            for y, m, v, sz in zip(c.ys, c.masses, c.volumes, c.sizes):
                b = next(bb for bb in self.bands if bb.lane == c.lane and bb.y == y)
                # a band split off by the ladder-guided fit sizes fine, but its volume is a model
                # estimate (faint shoulder of a bright band): not used as an ng standard
                if b.source != "ladder" and not sat[self.bands.index(b)] and m > 0 and v > 0:
                    xs.append(m * k)
                    ys.append(v)
                    ss.append(sz)
        if len(xs) < 3:
            return
        x, y, ls = np.asarray(xs, float), np.asarray(ys, float), np.log10(np.asarray(ss, float))
        slope = float((x * y).sum() / (x * x).sum())          # through origin
        n_sizes = len(np.unique(np.round(ls, 3)))
        model = self.params.mass_model if n_sizes >= 4 else "linear"
        deg = 3 if n_sizes >= 8 else 2 if n_sizes >= 5 else 1

        def fit(xx, yy, ll):
            if model == "size":
                cf = np.polyfit(ll, np.log10(yy / xx), deg)
                return lambda l_: 10 ** np.polyval(cf, np.clip(l_, ls.min(), ls.max())), cf
            sl = float((xx * yy).sum() / (xx * xx).sum())
            return lambda l_: np.full_like(np.asarray(l_, float), sl), None
        resp, coef = fit(x, y, ls)
        pred = x * resp(ls)
        ss_res = float(((y - pred) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum()) or 1.0
        loo = []
        for i in range(len(x)):
            keep = np.arange(len(x)) != i
            if keep.sum() >= deg + 2:
                r_i, _ = fit(x[keep], y[keep], ls[keep])
                loo.append(np.log10(y[i] / r_i(ls[i]) / x[i]))
        err = float(10 ** np.sqrt(np.mean(np.square(loo))) - 1) if loo else float("nan")
        unit = PRESETS.get(self.params.marker_preset or "", {}).get("mass_unit", "ng")
        self.mass_fit = {"model": model, "slope": slope, "coef": None if coef is None else [float(c) for c in coef],
                         "log_size_range": [float(ls.min()), float(ls.max())],
                         "load_scale": k, "r2": 1 - ss_res / ss_tot, "typical_error": err, "n": len(xs),
                         "min": float(x.min()), "max": float(x.max()), "unit": unit,
                         "points": [[float(a), float(b_), float(c)] for a, b_, c in zip(xs, ys, ss)]}

    def mass_for(self, volume: float, size: float) -> tuple[float, bool]:
        """Estimated mass of a band (needs :meth:`calibrate`). Returns (mass, size_corrected)."""
        mf = self.mass_fit
        if not mf:
            return float("nan"), False
        if mf["model"] == "size" and mf["coef"] and np.isfinite(size) and size > 0:
            lo, hi = mf["log_size_range"]
            r = 10 ** np.polyval(mf["coef"], np.clip(np.log10(size), lo, hi))
            return float(volume / r), True
        return float(volume / mf["slope"]), False

    def _size_for(self, lane: int, y: float) -> tuple[float, bool]:
        gel = (self.lanes[lane].tier, self.lanes[lane].column)
        # several gels in one image: only markers of the same gel describe this lane's migration
        cals = [c for c in self.calibrations if len(c.ys) >= 2
                and (self.lanes[c.lane].tier, self.lanes[c.lane].column) == gel]
        if not cals:
            return np.nan, False
        x = self.lanes[lane].x
        if len(cals) >= 2 and len({len(c.ys) for c in cals}) == 1 and \
                all(c.sizes == cals[0].sizes for c in cals):
            # two markers: interpolate each ladder band's position across the gel (corrects "smiling")
            cals = sorted(cals, key=lambda c: self.lanes[c.lane].x)
            xs = np.array([self.lanes[c.lane].x for c in cals])
            ys = np.array([c.ys for c in cals])
            yy = [float(np.interp(x, xs, ys[:, k])) if xs[0] <= x <= xs[-1] else
                  float(ys[0, k] + (ys[-1, k] - ys[0, k]) * (x - xs[0]) / (xs[-1] - xs[0]))
                  for k in range(ys.shape[1])]
            cal = MarkerCalibration(lane=lane, ys=yy, sizes=cals[0].sizes)
        else:
            cal = min(cals, key=lambda c: abs(self.lanes[c.lane].x - x))
        ls, ext = cal.log_size_at(y)
        return float(10 ** ls[0]), bool(ext[0])

    # --------------------------------------------------------- quantification
    def _volumes(self) -> list[float]:
        return self._volumes_info()[0]

    def _volumes_info(self) -> tuple[list[float], list[bool]]:
        """Background-corrected volume per band, and whether it came from peak fitting."""
        vols, fitted = [], []
        for b in self.bands:
            pr = self.lane_profile(b.lane)
            sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
            vols.append(float(pr["corr"][sel].sum()))
            fitted.append(False)
        if self.params.deconvolve:
            for group in self._overlap_groups():
                split = self._deconvolve(group)
                if split is not None:
                    for k, v in zip(group, split):
                        vols[k], fitted[k] = v, True
        return vols, fitted

    def _overlap_groups(self) -> list[list[int]]:
        """Runs of automatic bands in a lane whose integration ranges touch (split at a valley)."""
        groups, cur = [], []
        order = sorted(range(len(self.bands)), key=lambda k: (self.bands[k].lane, self.bands[k].y))
        for k in order:
            b = self.bands[k]
            prev = self.bands[cur[-1]] if cur else None
            touching = prev is not None and prev.lane == b.lane and b.y0 <= prev.y1 + 1.5
            if touching and b.source in AUTO_SOURCES and prev.source in AUTO_SOURCES:
                cur.append(k)
            else:
                if len(cur) > 1:
                    groups.append(cur)
                cur = [k] if b.source in AUTO_SOURCES else []
        if len(cur) > 1:
            groups.append(cur)
        return groups

    def _deconvolve(self, group: list[int]) -> list[float] | None:
        """Fit one Gaussian per band over the shared range; return the shared signal split by the
        fitted areas, or None if the fit is poor (then the valley split is kept)."""
        from scipy.optimize import curve_fit
        bands = [self.bands[k] for k in group]
        pr = self.lane_profile(bands[0].lane)
        y0, y1 = bands[0].y0, bands[-1].y1
        sel = (pr["y"] >= y0) & (pr["y"] <= y1)
        yy, cc = pr["y"][sel], pr["corr"][sel]
        n = len(bands)
        if len(yy) < 3 * n + 2 or cc.sum() <= 0:
            return None

        def model(y, *p):
            out = np.zeros_like(y, dtype=float)
            for j in range(n):
                a, mu, sg = p[3 * j: 3 * j + 3]
                out += a * np.exp(-0.5 * ((y - mu) / sg) ** 2)
            return out
        p0, lo, hi = [], [], []
        for b in bands:
            w = max(1.0, (b.y1 - b.y0) / 4)
            a0 = max(float(np.interp(b.y, yy, cc)), 1e-9)
            p0 += [a0, b.y, w]
            lo += [0.0, b.y - 3, 0.5]
            hi += [a0 * 3 + 1e-9, b.y + 3, max(2.0, (y1 - y0) / 2)]
        try:
            with np.errstate(all="ignore"):
                p, _ = curve_fit(model, yy, cc, p0=p0, bounds=(lo, hi), maxfev=4000)
        except (RuntimeError, ValueError):
            return None
        res = cc - model(yy, *p)
        r2 = 1 - float((res ** 2).sum()) / (float(((cc - cc.mean()) ** 2).sum()) or 1.0)
        if not np.isfinite(r2) or r2 < DECONV_MIN_R2:
            return None
        areas = np.array([p[3 * j] * p[3 * j + 2] for j in range(n)])
        if areas.sum() <= 0:
            return None
        return list(float(cc.sum()) * areas / areas.sum())

    def _saturated(self) -> list[int]:
        out = []
        for b in self.bands:
            pr = self.lane_profile(b.lane)
            sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
            out.append(int(pr["sat"][sel].sum()))
        return out

    def quantify(self) -> pd.DataFrame:
        if self.params.marker_lanes and self.params.marker_sizes:
            self.calibrate()
        vols, fitted = self._volumes_info()
        sats = self._saturated()
        marker_idx = {L - 1 for L in self.params.marker_lanes}
        rows = []
        lane_tot = {}
        for b, v in zip(self.bands, vols):
            lane_tot[b.lane] = lane_tot.get(b.lane, 0.0) + v
        cnt = {}
        def rf(b):
            lane = self.lanes[b.lane]
            well = lane.y0 + 10
            return (b.y - well) / max(1.0, lane.y1 - well)
        unit = self.params.marker_unit
        for b, v, s, fit in zip(self.bands, vols, sats, fitted):
            cnt[b.lane] = cnt.get(b.lane, 0) + 1
            pr = self.lane_profile(b.lane)
            sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
            size, ext = self._size_for(b.lane, b.y)
            r = {
                "lane": b.lane + 1, "lane_name": self.lanes[b.lane].name,
                "gel_row": self.lanes[b.lane].tier + 1,
                "gel_column": self.lanes[b.lane].column + 1,
                "is_marker": b.lane in marker_idx, "band": cnt[b.lane],
                "y_px": round(b.y, 1), "y0_px": round(b.y0, 1), "y1_px": round(b.y1, 1),
                "rf": round(rf(b), 4),
                f"size_{unit}": round(size, 1) if np.isfinite(size) else np.nan,
                "size_extrapolated": ext if np.isfinite(size) else None,
                "size_in_range": (bool(self.size_range()[0] <= size <= self.size_range()[1]) and not ext)
                if np.isfinite(size) else None,
                "volume": round(v, 1),
                "deconvolved": fit,
                "raw_volume": round(float(pr["raw"][sel].sum()), 1),
                "background": round(float(pr["baseline"][sel].sum()), 1),
                "peak_height": round(float(pr["corr"][sel].max()) if sel.any() else 0.0, 1),
                "pct_of_lane": round(100 * v / lane_tot[b.lane], 2) if lane_tot[b.lane] else np.nan,
                "saturated_px": s, "saturated": s > 0,
                "lane_coverage": round(self.coverage(b), 2), "source": b.source,
            }
            if self.mass_fit:
                m, corrected = self.mass_for(v, size)
                r[f"est_mass_{self.mass_fit['unit']}"] = round(m, 2)
                r["mass_size_corrected"] = corrected
                r["mass_out_of_range"] = not (self.mass_fit["min"] <= m <= self.mass_fit["max"])
                vol_ul = self.params.sample_loads_ul.get(str(b.lane + 1))
                if vol_ul:
                    r[f"conc_{self.mass_fit['unit']}_per_ul"] = round(m / float(vol_ul), 3)
            rows.append(r)
        df = pd.DataFrame(rows)
        df = self._model_columns(df)
        df = self._compare(df)
        n_sat = int(df["saturated"].sum()) if len(df) else 0
        extra = "; ".join(self.compare_warnings)
        self._log("quantify", {"norm_method": self.params.norm_method, "rna_integrity": self.params.rna_integrity},
                  f"{len(df)} bands quantified, {n_sat} saturated" + (f"; {extra}" if extra else ""))
        return df

    def _model_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """band_probability (band checker) and, after an AI review of these bands, ai_support."""
        if not len(df):
            return df
        from .checker import band_features, load_model, predict
        model = load_model()
        if model is not None:
            df["band_probability"] = np.round(predict(model, band_features(self)), 3)
        rev = getattr(self, "ai_result", None)
        if rev and rev.get("bands") == [(b.lane, round(b.y0), round(b.y1)) for b in self.bands]:
            df["ai_support"] = rev["support"]
        return df

    def ai_review(self, add_missing: bool = True) -> dict:
        """Second opinion of GelGenie's network (optional extra, see ai.py): flags autogel's bands
        it doesn't confirm and (optionally) adds the band regions it sees where autogel has none."""
        from . import ai
        ok, why = ai.available()
        if not ok:
            raise RuntimeError(f"AI second opinion not available: {why}")
        r = ai.review(self)
        added = 0
        if add_missing and r["missing"]:
            self.bands = sorted(self.bands + r["missing"], key=lambda b: (b.lane, b.y))
            added = len(r["missing"])
            r = ai.review(self)                      # support for the new list (added bands count as supported)
        self.ai_result = {"support": r["support"], "bands": [(b.lane, round(b.y0), round(b.y1)) for b in self.bands],
                          "model": ai.MODEL["name"]}
        n_un = sum(not s for s in r["support"])
        self._log("ai_review", {"model": ai.MODEL["name"], "add_missing": add_missing},
                  f"{len(self.bands)} bands: {n_un} not confirmed by the network, {added} added from it")
        return self.ai_result

    def _compare(self, df: pd.DataFrame) -> pd.DataFrame:
        """Per-channel volumes, band groups across lanes, normalisation, RNA integrity."""
        from . import compare
        self.compare_warnings: list[str] = []
        self.rna = pd.DataFrame()
        if not len(df):
            return df
        if len(self.gel.channels) > 1:              # every band measured in every channel
            for k, name in enumerate(self.gel.channel_names):
                col = compare.channel_column(name)
                vals = []
                for b in self.bands:
                    pr = self.lane_profile(b.lane, channel=k)
                    sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
                    vals.append(round(float(pr["corr"][sel].sum()), 1))
                df[col] = vals
        size_col = f"size_{self.params.marker_unit}"
        size_col = size_col if size_col in df and df.loc[~df["is_marker"], size_col].notna().any() else None
        df = compare.add_band_groups(df, size_col, self.params.group_tol)
        df, w = compare.normalise(df, self, size_col)
        self.compare_warnings += w
        if self.params.rna_integrity:
            self.rna = compare.rna_integrity(self, df, size_col)
        return df

    def lanes_table(self) -> pd.DataFrame:
        rows = []
        vols = self._volumes()
        for i, l in enumerate(self.lanes):
            pr = self.lane_profile(i)
            rows.append({"lane": i + 1, "name": l.name, "x_px": round(l.x, 1), "width_px": round(l.width, 1),
                         "y0_px": round(l.y0), "y1_px": round(l.y1),
                         "is_marker": (i + 1) in self.params.marker_lanes, "source": l.source,
                         "n_bands": sum(b.lane == i for b in self.bands),
                         "sum_band_volume": round(sum(v for b, v in zip(self.bands, vols) if b.lane == i), 1),
                         "total_lane_signal": round(float(pr["corr"].sum()), 1)})
        out = pd.DataFrame(rows)
        rna = getattr(self, "rna", None)
        if rna is not None and len(rna):
            out = out.merge(rna.drop(columns=["lane_name"]), on="lane", how="left")
        return out

    # ----------------------------------------------------------------- run all
    def run_auto(self) -> pd.DataFrame:
        self.detect_lanes()
        self.detect_bands()
        return self.quantify()

    # --------------------------------------------------------------- session
    def to_dict(self) -> dict:
        return {"autogel_version": __version__, "image": self.path,
                "orientation": self.gel.orientation, "params": asdict(self.params),
                "well_y": self.well_y, "ladder_refined": sorted(self._ladder_refined),
                "lanes": [asdict(l) for l in self.lanes],
                "bands": [asdict(b) for b in self.bands], "log": self.log}

    def save_session(self, path: str):
        def conv(o):
            if isinstance(o, (np.floating, np.integer)):
                return o.item()
            if isinstance(o, np.bool_):
                return bool(o)
            raise TypeError(type(o))
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=1, default=conv)

    @classmethod
    def load_session(cls, path: str, image_path: str | None = None) -> GelAnalysis:
        """Rebuild an analysis exactly (same orientation, lanes and bands, incl. manual edits)."""
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        params = Params(**d["params"])
        ga = cls(image_path or d["image"], params, orientation=d["orientation"])
        ga.log = d.get("log", []) + ga.log[-1:]
        ga.well_y = d.get("well_y", 0.0)
        ga.lanes = [Lane(**l) for l in d["lanes"]]          # older sessions: no "bend" -> straight
        ga.bands = [Band(**b) for b in d["bands"]]
        # sessions before 0.5: treat marker lanes as refined (the user has seen and edited them)
        ga._ladder_refined = set(d.get("ladder_refined", [L - 1 for L in params.marker_lanes]))
        ga._log("load_session", {"session": os.path.abspath(path)},
                f"{len(ga.lanes)} lanes, {len(ga.bands)} bands restored")
        return ga

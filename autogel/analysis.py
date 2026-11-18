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
from .imageio import GelImage, load_gel
from .markers import PRESETS, MarkerCalibration, assign_marker


# ------------------------------------------------------------------ data model
@dataclass
class Params:
    invert: bool | None = None          # None = auto (bright bands on dark = False)
    n_lanes: int | None = None          # None = auto (per gel row)
    n_tiers: int = 1                    # gel rows = combs stacked in one image
    tier_rows: list[float] | None = None  # split rows between gel rows (None = automatic)
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


@dataclass
class Lane:
    x: float
    width: float
    y0: float
    y1: float
    name: str = ""
    source: str = "auto"
    tier: int = 0                       # gel row (0 = top comb)

    @property
    def x0(self):
        return self.x - self.width / 2

    @property
    def x1(self):
        return self.x + self.width / 2


@dataclass
class Band:
    lane: int        # 0-based index into lanes
    y: float         # peak row
    y0: float        # integration start row
    y1: float        # integration end row (inclusive)
    source: str = "auto"


class GelAnalysis:
    """Holds the full state of one gel analysis plus a log of every step."""

    def __init__(self, path: str, params: Params | None = None, orientation: dict | None = None):
        self.params = params or Params()
        self.log: list[dict] = []
        self.listeners: list = []
        self.gel: GelImage = load_gel(path, orientation=orientation, invert=self.params.invert)
        self.path = os.path.abspath(path)
        self.lanes: list[Lane] = []
        self.bands: list[Band] = []
        self.well_y: float = 0.0
        self.calibrations: list[MarkerCalibration] = []
        self.mass_fit: dict | None = None
        self._profiles: dict[int, dict] = {}
        H = self.gel.data.shape[0]
        if self.params.bg_radius is None:
            self.params.bg_radius = max(15.0, H / 40)
        self._log("load", {"path": self.path, "orientation": self.gel.orientation,
                           "dtype": self.gel.dtype, "shape": list(self.gel.data.shape),
                           "saturation_level": self.gel.saturation_level},
                  "; ".join(self.gel.warnings) or "ok")

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

    # ------------------------------------------------------------------ lanes
    def detect_lanes(self) -> list[Lane]:
        """Find lanes. With several gel rows (stacked gels / two combs) every row is handled
        separately: its own wells, lane comb and lane names (A1.., B1..)."""
        a = self.data
        H, W = a.shape
        # Band-like features only: a vertical top-hat keeps short horizontal
        # structures (bands, wells) and removes long vertical ones (tray edges,
        # smears, illumination gradients).
        bandimg = self._band_feature_image()
        tiers = self.tier_bounds()
        lanes: list[Lane] = []
        info = []
        for t, (ra, rb) in enumerate(tiers):
            h = rb - ra
            prof = bandimg[ra + int(h * .02): rb - int(h * .02)].mean(0)
            # baseline window must be wider than any lane (wide lanes were split in two otherwise)
            prof = prof - ndimage.minimum_filter1d(ndimage.gaussian_filter1d(prof, 2), max(5, W // 4))
            comb = prof / (prof.max() + 1e-12)
            peaks, spacing = self._fit_lane_grid(comb, self.params.n_lanes)
            width = 0.8 * spacing
            well = ra + self._find_well(a[ra:rb], peaks, width)
            y0 = max(float(ra), well - 10)
            y1 = float(rb - 1 - (int(H * .005) if t == len(tiers) - 1 else 2))
            prefix = chr(65 + t) if len(tiers) > 1 else ""
            lanes += [Lane(x=float(p), width=float(width), y0=y0, y1=y1, name=f"{prefix}{i + 1}", tier=t)
                      for i, p in enumerate(peaks)]
            info.append(f"{prefix + ': ' if prefix else ''}{len(peaks)} lanes, width {width:.0f} px, "
                        f"spacing {spacing:.1f} px, well row {well:.0f}")
        self.lanes = lanes
        self.well_y = lanes[0].y0 + 10 if lanes else 0.0
        self.bands = []
        self._profiles = {}
        self._log("detect_lanes", {"n_lanes_param": self.params.n_lanes, "n_tiers": len(tiers),
                                   "tier_bounds": [[int(x), int(y)] for x, y in tiers]},
                  "; ".join(info))
        return self.lanes

    def _find_well(self, a: np.ndarray, peaks, width) -> float:
        """Well row inside one gel row: every lane has a well, only some have bands -> median
        over lanes, then the TOPMOST clear peak in the top 20 % (not the strongest one)."""
        H = a.shape[0]
        top = a[: max(5, int(H * .2))]
        per_lane = [top[:, max(0, int(p - width / 2)): int(p + width / 2) + 1].mean(1) for p in peaks]
        rowp = np.median(per_lane, axis=0) if per_lane else top.mean(1)
        rowp = ndimage.gaussian_filter1d(rowp, 2)
        pk, pr = signal.find_peaks(rowp, prominence=0)
        noise = 3 * np.std(np.diff(rowp))
        good = []
        if len(pk):
            thr = max(noise, 0.15 * pr["prominences"].max())
            good = [k for k, pp in zip(pk, pr["prominences"]) if pp > thr]
        return float(good[0]) if good else float(int(H * .01))

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
                x = old_lanes[b.lane].x
                hit = [j for j, l in enumerate(lanes) if l.x0 <= x <= l.x1]
                if hit:
                    b.lane = hit[0]
                    new_bands.append(b)
        self.bands = new_bands
        self._log("set_lanes", {"reason": reason,
                                "lanes": [[round(l.x, 1), round(l.width, 1), round(l.y0), round(l.y1)] for l in lanes]},
                  f"{old} -> {len(lanes)} lanes, {len(new_bands)} bands kept")

    def set_lane_names(self, names: list[str]):
        for l, n in zip(self.lanes, names):
            l.name = n
        self._log("set_lane_names", {"names": names})

    # ------------------------------------------------------------- profiles
    def lane_profile(self, i: int) -> dict:
        if i in self._profiles:
            return self._profiles[i]
        lane = self.lanes[i]
        a = self.data
        H, W = a.shape
        hw = lane.width * self.params.lane_width_frac / 2
        c0, c1 = max(0, int(round(lane.x - hw))), min(W, int(round(lane.x + hw)) + 1)
        ya, yb = sorted((lane.y0, lane.y1))
        r0, r1 = max(0, int(ya)), min(H, int(yb) + 1)
        r1 = max(r0, r1)
        strip = a[r0:r1, c0:c1]
        raw = strip.sum(1, dtype=np.float64)
        sat = self.gel.raw_max[r0:r1, c0:c1].sum(1)
        base = self._baseline(raw, i, (r0, r1, c0, c1))
        prof = {"y": np.arange(r0, r1, dtype=float), "raw": raw, "baseline": base,
                "corr": raw - base, "sat": sat, "cols": (c0, c1)}
        self._profiles[i] = prof
        return prof

    def _baseline(self, raw: np.ndarray, i: int, box) -> np.ndarray:
        m = self.params.bg_method
        if m == "none":
            return np.zeros_like(raw)
        if m == "lane_edges":
            r0, r1, c0, c1 = box
            a = self.data
            g = max(2, int(self.lanes[i].width * 0.1))
            left = a[r0:r1, max(0, c0 - 2 * g): max(0, c0 - g)]
            right = a[r0:r1, c1 + g: c1 + 2 * g]
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
        weak = 0
        if new and self.params.min_gel_rel > 0:
            others = [self._peak_height(b) for b in self.bands]
            ref = float(np.percentile([h for _, h in new] + others, 95))
            keep = [(b, h) for b, h in new if h >= self.params.min_gel_rel * ref]
            weak = len(new) - len(keep)
            new = keep
        for b, _ in new:
            self.bands.append(b)
            total += 1
        self.bands.sort(key=lambda b: (b.lane, b.y))
        self._log("detect_bands", {"lanes": None if lanes is None else [l + 1 for l in lanes],
                                   "bg_method": self.params.bg_method, "bg_radius": self.params.bg_radius,
                                   "min_snr": self.params.min_snr, "min_rel_height": self.params.min_rel_height,
                                   "min_gel_rel": self.params.min_gel_rel,
                                   "skip_well_px": self.params.skip_well_px},
                  f"{total} bands ({rejected} rejected as narrower than "
                  f"{self.params.min_coverage:.0%} of lane width: dust/specks"
                  + (f"; {weak} below {self.params.min_gel_rel:.0%} of the gel's strong bands" if weak else "")
                  + ")")
        return self.bands

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
        y0, y1 = int(b.y0), int(b.y1) + 1
        block = self.data[y0:y1, c0:c1].sum(0, dtype=np.float64)
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

    def calibrate(self) -> list[MarkerCalibration]:
        p = self.params
        self.calibrations, self.mass_fit = [], None
        if not p.marker_sizes or not p.marker_lanes:
            return []
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
                  (f"; mass fit slope {self.mass_fit['slope']:.4g} vol/{self.mass_fit['unit']}, "
                   f"R²={self.mass_fit['r2']:.3f}, n={self.mass_fit['n']}" if self.mass_fit else ""))
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
        xs, ys = [], []
        k = self.mass_scale()
        sat = self._saturated()
        for c in self.calibrations:
            if not c.masses:
                continue
            for y, m, v in zip(c.ys, c.masses, c.volumes):
                b = next(bb for bb in self.bands if bb.lane == c.lane and bb.y == y)
                if not sat[self.bands.index(b)]:
                    xs.append(m * k)
                    ys.append(v)
        if len(xs) < 3:
            return
        x, y = np.asarray(xs, float), np.asarray(ys, float)
        slope = float((x * y).sum() / (x * x).sum())          # through origin
        pred = slope * x
        ss_res = float(((y - pred) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum()) or 1.0
        unit = PRESETS.get(self.params.marker_preset or "", {}).get("mass_unit", "ng")
        self.mass_fit = {"slope": slope, "load_scale": k, "r2": 1 - ss_res / ss_tot, "n": len(xs),
                         "min": float(x.min()), "max": float(x.max()), "unit": unit,
                         "points": list(zip(xs, ys))}

    def _size_for(self, lane: int, y: float) -> tuple[float, bool]:
        tier = self.lanes[lane].tier
        # stacked gels: only markers of the same gel row describe this lane's migration
        cals = [c for c in self.calibrations if len(c.ys) >= 2 and self.lanes[c.lane].tier == tier]
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
        out = []
        for b in self.bands:
            pr = self.lane_profile(b.lane)
            sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
            out.append(float(pr["corr"][sel].sum()))
        return out

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
        vols = self._volumes()
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
        for b, v, s in zip(self.bands, vols, sats):
            cnt[b.lane] = cnt.get(b.lane, 0) + 1
            pr = self.lane_profile(b.lane)
            sel = (pr["y"] >= b.y0) & (pr["y"] <= b.y1)
            size, ext = self._size_for(b.lane, b.y)
            r = {
                "lane": b.lane + 1, "lane_name": self.lanes[b.lane].name,
                "gel_row": self.lanes[b.lane].tier + 1,
                "is_marker": b.lane in marker_idx, "band": cnt[b.lane],
                "y_px": round(b.y, 1), "y0_px": round(b.y0, 1), "y1_px": round(b.y1, 1),
                "rf": round(rf(b), 4),
                f"size_{unit}": round(size, 1) if np.isfinite(size) else np.nan,
                "size_extrapolated": ext if np.isfinite(size) else None,
                "size_in_range": (bool(self.size_range()[0] <= size <= self.size_range()[1]) and not ext)
                if np.isfinite(size) else None,
                "volume": round(v, 1),
                "raw_volume": round(float(pr["raw"][sel].sum()), 1),
                "background": round(float(pr["baseline"][sel].sum()), 1),
                "peak_height": round(float(pr["corr"][sel].max()) if sel.any() else 0.0, 1),
                "pct_of_lane": round(100 * v / lane_tot[b.lane], 2) if lane_tot[b.lane] else np.nan,
                "saturated_px": s, "saturated": s > 0,
                "lane_coverage": round(self.coverage(b), 2), "source": b.source,
            }
            if self.mass_fit:
                m = v / self.mass_fit["slope"]
                r[f"est_mass_{self.mass_fit['unit']}"] = round(m, 2)
                r["mass_out_of_range"] = not (self.mass_fit["min"] <= m <= self.mass_fit["max"])
                vol_ul = self.params.sample_loads_ul.get(str(b.lane + 1))
                if vol_ul:
                    r[f"conc_{self.mass_fit['unit']}_per_ul"] = round(m / float(vol_ul), 3)
            rows.append(r)
        df = pd.DataFrame(rows)
        n_sat = int(df["saturated"].sum()) if len(df) else 0
        self._log("quantify", {}, f"{len(df)} bands quantified, {n_sat} saturated")
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
        return pd.DataFrame(rows)

    # ----------------------------------------------------------------- run all
    def run_auto(self) -> pd.DataFrame:
        self.detect_lanes()
        self.detect_bands()
        return self.quantify()

    # --------------------------------------------------------------- session
    def to_dict(self) -> dict:
        return {"autogel_version": __version__, "image": self.path,
                "orientation": self.gel.orientation, "params": asdict(self.params),
                "well_y": self.well_y,
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
        ga.lanes = [Lane(**l) for l in d["lanes"]]
        ga.bands = [Band(**b) for b in d["bands"]]
        ga._log("load_session", {"session": os.path.abspath(path)},
                f"{len(ga.lanes)} lanes, {len(ga.bands)} bands restored")
        return ga

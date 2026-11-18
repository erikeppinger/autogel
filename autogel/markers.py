"""Marker (ladder) presets and size / mass calibration.

Sizes are listed from the top of the gel (largest) to the bottom (smallest).
``masses`` (optional) are the ng (DNA) / µg (protein) per band for the
standard load given in the vendor datasheet and enable a mass standard curve.

Ladders live in ``ladders.json`` (built-in) and ``~/.autogel/ladders.json``
(user-defined, overrides built-in entries of the same name).
ALWAYS check the presets against the datasheet of the lot you use.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator


def _user_library_path() -> str:
    base = os.environ.get("AUTOGEL_HOME") or os.path.join(os.path.expanduser("~"), ".autogel")
    return os.path.join(base, "ladders.json")


def load_library() -> dict[str, dict]:
    """Ladder library, in this order (later files override earlier ones by name):

    1. ``ladders.json``      curated, checked against vendor datasheets (with masses)
    2. ``ladders_ape.json``  imported from ApE (sizes only, unverified); skipped when the same
                             catalog number is already in (1)
    3. ``~/.autogel/ladders.json``  the user's own ladders
    """
    here = os.path.dirname(__file__)
    out: dict[str, dict] = {}
    curated_catalogs: set[str] = set()
    sources = [(os.path.join(here, "ladders.json"), "builtin"),
               (os.path.join(here, "ladders_ape.json"), "ape"),
               (_user_library_path(), "user")]
    for path, kind in sources:
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        for e in data.get("ladders", []):
            e = dict(e)
            if kind == "ape" and e.get("catalog") in curated_catalogs:
                continue
            if kind == "builtin" and e.get("catalog"):
                curated_catalogs.add(e["catalog"])
            e["user_defined"] = kind == "user"
            if len(e.get("masses") or []) not in (0, len(e["sizes"])):
                e["masses"] = None
            out[e["name"]] = e
    return out


def save_user_ladder(name: str, sizes: list[float], unit: str = "bp",
                     masses: list[float] | None = None, source: str = "", note: str = "") -> str:
    """Add or replace a ladder in the user library; returns the file path."""
    if masses is not None and len(masses) != len(sizes):
        raise ValueError(f"{len(sizes)} sizes but {len(masses)} masses")
    path = _user_library_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = {"format": 1, "ladders": []}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    data["ladders"] = [e for e in data.get("ladders", []) if e.get("name") != name]
    data["ladders"].append({"name": name, "unit": unit, "sizes": sorted(sizes, reverse=True),
                            "masses": masses, "source": source, "note": note, "verified": True})
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    PRESETS.clear()
    PRESETS.update(load_library())
    return path


PRESETS: dict[str, dict] = load_library()
# Approximate effective separation ranges (common lab guidelines; adjust to your gel system,
# buffer and run length). Used only as suggestions for the "valid size range".
RESOLUTION_GUIDE: dict[str, tuple[float, float, str]] = {
    "agarose 0.5 %": (1000, 30000, "bp"),
    "agarose 0.7 %": (800, 12000, "bp"),
    "agarose 1.0 %": (500, 10000, "bp"),
    "agarose 1.2 %": (400, 7000, "bp"),
    "agarose 1.5 %": (200, 3000, "bp"),
    "agarose 2.0 %": (50, 2000, "bp"),
    "SDS-PAGE 7.5 %": (40, 250, "kDa"),
    "SDS-PAGE 10 %": (20, 150, "kDa"),
    "SDS-PAGE 12 %": (15, 100, "kDa"),
    "SDS-PAGE 15 %": (10, 70, "kDa"),
    "SDS-PAGE 4-20 % gradient": (10, 250, "kDa"),
}

DEFAULT_PRESET = "NEB Quick-Load Purple 1 kb Plus (N0550)"


def parse_sizes(text: str) -> list[float]:
    """Parse a user list like '10000, 8000; 6000 5000' into floats (sorted desc)."""
    vals = [float(t) for t in text.replace(";", ",").replace(" ", ",").split(",") if t.strip()]
    return sorted(vals, reverse=True)


@dataclass
class MarkerCalibration:
    """Maps migration position (row) -> size for one marker lane."""
    lane: int
    ys: list[float]            # band rows, top -> bottom
    sizes: list[float]         # matching sizes, large -> small
    masses: list[float] | None = None
    volumes: list[float] | None = None  # measured bg-corrected volumes
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self):
        order = np.argsort(self.ys)
        self.ys = [float(self.ys[i]) for i in order]
        self.sizes = [float(self.sizes[i]) for i in order]
        if self.masses is not None:
            self.masses = [float(self.masses[i]) for i in order]
        if self.volumes is not None:
            self.volumes = [float(self.volumes[i]) for i in order]
        logs = np.log10(self.sizes)
        if np.any(np.diff(logs) >= 0):
            self.warnings.append(
                f"Marker lane {self.lane}: sizes are not strictly decreasing with migration "
                "- check the band assignment.")

    def log_size_at(self, y: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        """Return (log10 size, extrapolated flag). PCHIP inside, linear outside."""
        y = np.atleast_1d(np.asarray(y, float))
        ys = np.asarray(self.ys)
        ls = np.log10(self.sizes)
        out = np.empty_like(y)
        inside = (y >= ys[0]) & (y <= ys[-1])
        if len(ys) >= 2:
            out[inside] = PchipInterpolator(ys, ls)(y[inside])
            lo = y < ys[0]
            hi = y > ys[-1]
            s_lo = (ls[1] - ls[0]) / (ys[1] - ys[0])
            s_hi = (ls[-1] - ls[-2]) / (ys[-1] - ys[-2])
            out[lo] = ls[0] + s_lo * (y[lo] - ys[0])
            out[hi] = ls[-1] + s_hi * (y[hi] - ys[-1])
        else:
            out[:] = np.nan
        return out, ~inside

    def to_dict(self):
        return {"lane": self.lane, "ys": self.ys, "sizes": self.sizes,
                "masses": self.masses, "volumes": self.volumes, "warnings": self.warnings}


def _fit_quality(ys, sizes, vols=None, masses=None) -> float:
    """Lower is better: roughness of log10(size) vs migration (cubic fit residual, relative to
    the ladder's span) plus, if masses are known, how badly band brightness disagrees with them."""
    y = np.asarray(ys, float)
    ls = np.log10(np.asarray(sizes, float))
    if len(y) < 5:
        return 0.0
    yn = (y - y.mean()) / (y.std() + 1e-9)
    coef = np.polyfit(yn, ls, 3)
    rms = float(np.sqrt(np.mean((np.polyval(coef, yn) - ls) ** 2)))
    score = rms / (ls.max() - ls.min() + 1e-9)
    if np.any(np.diff(ls) >= 0):
        score += 10.0
    if vols is not None and masses is not None and len(vols) == len(masses) >= 4:
        r = np.corrcoef(np.asarray(vols, float), np.asarray(masses, float))[0, 1]
        if np.isfinite(r):
            score += 0.05 * (1 - r)
    return score


def assign_marker(band_ys: list[float], sizes: list[float], masses: list[float] | None = None,
                  volumes: list[float] | None = None, max_skip: int = 3
                  ) -> tuple[list[float], list[float], list[float] | None, list[str]]:
    """Pair detected marker bands with ladder sizes.

    Equal counts -> direct pairing (top -> bottom). If counts differ by up to ``max_skip``,
    every way of dropping the surplus detected bands (dust, artefacts) or ladder sizes (bands
    not resolved / run off) is tried, and the assignment giving the smoothest size-vs-migration
    curve (and, with known masses, the best brightness match) wins. A warning always says
    what was done so the user can check it.
    """
    order = np.argsort(band_ys)
    ys = [float(band_ys[i]) for i in order]
    vols = [float(volumes[i]) for i in order] if volumes is not None else None
    sizes = list(sizes)
    warnings: list[str] = []
    m, n = len(ys), len(sizes)
    if m == n:
        return ys, sizes, (list(masses) if masses else None), warnings
    from itertools import combinations
    k = abs(m - n)
    if k > max_skip or min(m, n) < 3:
        q = min(m, n)
        warnings.append(f"Marker: {m} bands detected but the ladder has {n} sizes. Paired the top {q} "
                        "only - please correct the marker bands.")
        return ys[:q], sizes[:q], (list(masses[:q]) if masses else None), warnings
    best = None
    if m > n:     # drop k detected bands
        for drop in combinations(range(m), k):
            keep = [i for i in range(m) if i not in drop]
            sc = _fit_quality([ys[i] for i in keep], sizes, [vols[i] for i in keep] if vols else None, masses)
            if best is None or sc < best[0]:
                best = (sc, keep, list(range(n)))
        _, keep_y, keep_s = best
        dropped = [ys[i] for i in range(m) if i not in keep_y]
        warnings.append(f"Marker: {m} bands detected, ladder has {n}. Ignored the band(s) at row(s) "
                        f"{', '.join(f'{y:.0f}' for y in dropped)} (best fit) - please check.")
    else:         # drop k ladder sizes (missing bands)
        for drop in combinations(range(n), k):
            keep = [i for i in range(n) if i not in drop]
            sc = _fit_quality(ys, [sizes[i] for i in keep], vols,
                              [masses[i] for i in keep] if masses else None)
            if best is None or sc < best[0]:
                best = (sc, list(range(m)), keep)
        _, keep_y, keep_s = best
        missing = [sizes[i] for i in range(n) if i not in keep_s]
        warnings.append(f"Marker: only {m} of {n} ladder bands found; assumed missing: "
                        f"{', '.join(f'{s:g}' for s in missing)} (best fit) - please check.")
    return ([ys[i] for i in keep_y], [sizes[i] for i in keep_s],
            [masses[i] for i in keep_s] if masses else None, warnings)

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


def assign_marker(band_ys: list[float], sizes: list[float],
                  masses: list[float] | None = None) -> tuple[list[float], list[float], list[float] | None, list[str]]:
    """Pair detected marker bands with ladder sizes.

    Equal counts -> direct pairing (top->bottom). Otherwise the bands are
    paired from the top and a warning is returned; the user should then
    delete/add marker bands in the GUI so that the counts match.
    """
    ys = sorted(band_ys)
    warnings = []
    n = min(len(ys), len(sizes))
    if len(ys) != len(sizes):
        warnings.append(
            f"Marker: {len(ys)} bands detected but the ladder has {len(sizes)} sizes. "
            f"Paired the top {n} only - please correct the marker bands.")
    m = masses[:n] if masses is not None else None
    return ys[:n], list(sizes[:n]), m, warnings

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
import math
import os
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator

REPO = "erikeppinger/autogel"
LIBRARY_URL = f"https://raw.githubusercontent.com/{REPO}/main/autogel/ladders.json"
ISSUE_URL = f"https://github.com/{REPO}/issues/new"
MAX_LIBRARY_BYTES = 5_000_000
_INTERNAL = ("user_defined", "library")
UNITS = ("bp", "nt", "kDa")              # DNA base pairs, RNA nucleotides, protein kDa


def _home() -> str:
    return os.environ.get("AUTOGEL_HOME") or os.path.join(os.path.expanduser("~"), ".autogel")


def _user_library_path() -> str:
    return os.path.join(_home(), "ladders.json")


def _update_library_path() -> str:
    return os.path.join(_home(), "ladders_update.json")


def _read(path: str) -> dict:
    if not os.path.exists(path):
        return {"format": 2, "ladders": []}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _write(path: str, data: dict):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, path)          # never leave a half-written library behind


def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and np.isfinite(x)


def validate_ladder(e: dict) -> list[str]:
    """Problems with one ladder entry (empty list = fine). Used for import, update and saving."""
    if not isinstance(e, dict):
        return ["entry is not an object"]
    err = []
    name = e.get("name")
    if not isinstance(name, str) or not name.strip():
        err.append("missing name")
    if e.get("unit") not in UNITS:
        err.append(f"unit must be bp, nt or kDa, not {e.get('unit')!r}")
    sizes = e.get("sizes")
    if not isinstance(sizes, list) or len(sizes) < 2 or not all(_is_num(x) and x > 0 for x in sizes):
        err.append("sizes must be a list of at least 2 positive numbers")
        sizes = None
    elif any(b >= a for a, b in zip(sizes, sizes[1:])):
        err.append("sizes must be listed top to bottom (strictly decreasing, no duplicates)")
    masses = e.get("masses")
    if masses not in (None, []):
        if not isinstance(masses, list) or not all(_is_num(x) and x >= 0 for x in masses):
            err.append("masses must be non-negative numbers")
        elif sizes is not None and len(masses) != len(sizes):
            err.append(f"{len(sizes)} sizes but {len(masses)} masses")
    ref = e.get("reference_load")
    if ref is not None and not (isinstance(ref, dict) and _is_num(ref.get("amount")) and ref["amount"] > 0):
        err.append("reference_load must be like {'amount': 1.0, 'unit': 'µg'}")
    for k in ("source", "vendor", "catalog", "note"):
        if e.get(k) is not None and not isinstance(e.get(k), str):
            err.append(f"{k} must be text")
    return err


def load_library() -> dict[str, dict]:
    """Ladder library, in this order (later files override earlier ones by name):

    1. ``ladders.json``      curated, checked against vendor datasheets (with masses)
    2. ``~/.autogel/ladders_update.json``  newer curated library, downloaded on request with
                             :func:`update_library` (between releases)
    3. ``ladders_ape.json``  imported from ApE (sizes only, unverified); skipped when the same
                             catalog number is already in (1) or (2)
    4. ``~/.autogel/ladders.json``  the user's own ladders
    Invalid entries are skipped, so one bad entry can't break the program.
    """
    here = os.path.dirname(__file__)
    out: dict[str, dict] = {}
    curated_catalogs: set[str] = set()
    sources = [(os.path.join(here, "ladders.json"), "builtin"),
               (_update_library_path(), "update"),
               (os.path.join(here, "ladders_ape.json"), "ape"),
               (_user_library_path(), "user")]
    for path, kind in sources:
        if not os.path.exists(path):
            continue
        try:
            data = _read(path)
        except (OSError, ValueError):
            continue
        for e in data.get("ladders", []) if isinstance(data, dict) else []:
            if validate_ladder(e):
                continue
            e = dict(e)
            if kind == "ape" and e.get("catalog") in curated_catalogs:
                continue
            if kind in ("builtin", "update") and e.get("catalog"):
                curated_catalogs.add(e["catalog"])
            e["user_defined"] = kind == "user"
            e["library"] = kind
            if not e.get("masses"):
                e["masses"] = None
            out[e["name"]] = e
    return out


def _reload():
    PRESETS.clear()
    PRESETS.update(load_library())


def save_user_ladder(name: str, sizes: list[float], unit: str = "bp",
                     masses: list[float] | None = None, source: str = "", note: str = "",
                     vendor: str = "", catalog: str = "", reference_load: dict | None = None) -> str:
    """Add or replace a ladder in the user library; returns the file path.

    Sizes may be typed in any order; masses stay with their size."""
    if masses and len(masses) != len(sizes):
        raise ValueError(f"{len(sizes)} sizes but {len(masses)} masses")
    pairs = sorted(zip(sizes, masses or [None] * len(sizes)), key=lambda t: -float(t[0]))
    e = {"name": (name or "").strip(), "vendor": vendor, "catalog": catalog, "unit": unit,
         "sizes": [float(s) for s, _ in pairs],
         "masses": [float(m) for _, m in pairs] if masses else None,
         "mass_unit": "µg" if unit == "kDa" else "ng",
         "reference_load": reference_load, "source": source, "note": note, "verified": False}
    errs = validate_ladder(e)
    if errs:
        raise ValueError("; ".join(errs))
    path = _user_library_path()
    data = _read(path)
    data["ladders"] = [x for x in data.get("ladders", []) if x.get("name") != e["name"]] + [e]
    _write(path, data)
    _reload()
    return path


def delete_user_ladder(name: str) -> bool:
    """Remove one of the user's own ladders (built-in ones can't be deleted)."""
    path = _user_library_path()
    data = _read(path)
    keep = [x for x in data.get("ladders", []) if x.get("name") != name]
    if len(keep) == len(data.get("ladders", [])):
        return False
    data["ladders"] = keep
    _write(path, data)
    _reload()
    return True


def export_ladders(names: list[str], path: str) -> int:
    """Write the named ladders to a file another lab can import (or attach to a GitHub issue)."""
    lad = [{k: v for k, v in PRESETS[n].items() if k not in _INTERNAL} for n in names if n in PRESETS]
    _write(path, {"format": 2, "ladders": lad})
    return len(lad)


def import_ladders(path: str, overwrite: bool = False) -> dict:
    """Add the ladders of a file to the user library. Returns {added, replaced, skipped, errors}.

    Accepts a library file ({"ladders": [...]}), a list of ladders or a single ladder object."""
    res: dict[str, list] = {"added": [], "replaced": [], "skipped": [], "errors": []}
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError) as e:
        res["errors"].append(f"cannot read {os.path.basename(path)}: {e}")
        return res
    if isinstance(data, dict) and "sizes" in data:
        entries = [data]
    elif isinstance(data, dict):
        entries = data.get("ladders", [])
    else:
        entries = data
    if not isinstance(entries, list) or not entries:
        res["errors"].append("no ladders found in the file")
        return res
    upath = _user_library_path()
    lib = _read(upath)
    mine = {x.get("name"): x for x in lib.get("ladders", [])}
    for e in entries:
        errs = validate_ladder(e)
        if errs:
            res["errors"].append(f"{e.get('name', '?') if isinstance(e, dict) else '?'}: {'; '.join(errs)}")
            continue
        e = {k: v for k, v in e.items() if k not in _INTERNAL}
        e["verified"] = False                               # imported = not checked by us
        if e["name"] in mine or e["name"] in PRESETS:
            if not overwrite:
                res["skipped"].append(e["name"])
                continue
            res["replaced"].append(e["name"])
        else:
            res["added"].append(e["name"])
        mine[e["name"]] = e
    if res["added"] or res["replaced"]:
        lib["ladders"] = list(mine.values())
        _write(upath, lib)
        _reload()
    return res


def update_library(url: str = LIBRARY_URL, timeout: float = 20.0) -> tuple[int, str]:
    """Download the current curated library from the project's GitHub page (only when the user
    asks for it). Every entry is validated; an invalid or oversized file is rejected and the
    old library kept. Returns (number of ladders, saved path)."""
    import datetime as _dt
    import urllib.request
    if not url.startswith("https://"):
        raise ValueError("library URL must start with https://")
    import ssl
    try:        # python.org's macOS Python ships without root certificates; certifi has them
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        ctx = ssl.create_default_context()
    req = urllib.request.Request(url, headers={"User-Agent": "autogel"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:  # noqa: S310 - https only
        raw = r.read(MAX_LIBRARY_BYTES + 1)
    return _store_update(raw, url, _dt.datetime.now().isoformat(timespec="seconds"))


def _store_update(raw: bytes, url: str, when: str) -> tuple[int, str]:
    if len(raw) > MAX_LIBRARY_BYTES:
        raise ValueError("downloaded library is too large")
    data = json.loads(raw.decode("utf-8"))
    entries = data.get("ladders") if isinstance(data, dict) else None
    if not isinstance(entries, list) or not entries:
        raise ValueError("downloaded file contains no ladders")
    bad = [f"{e.get('name', '?') if isinstance(e, dict) else '?'}: {'; '.join(v)}"
           for e in entries if (v := validate_ladder(e))]
    if bad:
        raise ValueError("downloaded library has invalid entries: " + " | ".join(bad[:3]))
    data["fetched"], data["fetched_from"] = when, url
    path = _update_library_path()
    _write(path, data)
    _reload()
    return len(entries), path


def contribution_url(e: dict) -> str:
    """Pre-filled 'suggest a ladder' issue on GitHub, opened in the browser by the user; nothing
    is sent automatically. Maintainers check it against the datasheet and add it to ladders.json,
    so it ships with the next release and reaches existing installs via 'Update library'."""
    from urllib.parse import urlencode

    def fmt(v):
        return ", ".join(f"{x:g}" for x in v) if v else ""
    ref = e.get("reference_load") or {}
    entry = {k: v for k, v in e.items() if k not in _INTERNAL}
    q = {"template": "ladder.yml", "title": f"Ladder: {e.get('name', '')}",
         "name": e.get("name", ""), "vendor": e.get("vendor", ""), "catalog": e.get("catalog", ""),
         "unit": e.get("unit", "bp"), "sizes": fmt(e.get("sizes")), "masses": fmt(e.get("masses")),
         "reference_load": f"{ref.get('amount', '')} {ref.get('unit', '')}".strip(),
         "source": e.get("source", ""), "json": json.dumps(entry, ensure_ascii=False)}
    return ISSUE_URL + "?" + urlencode(q)


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


# Ladder alignment (D29). Two neighbouring steps of a real ladder have nearly the same spacing
# per log(size) (pixels per decade change by ~20 % from step to step; 95 % within 50 %) and
# nearly the same brightness per ng (~15 %). A merged pair, a missed band, a well or dye-front
# band breaks one of these locally - also where a global curve would just bend around it.
ALIGN_RATE_SD = 0.4        # log ratio of consecutive spacings counted as one "standard deviation"
ALIGN_BRIGHT_SD = 0.2      # same for log(volume / ng) between neighbours
ALIGN_BRIGHT_W = 0.25      # weight of the brightness term (when masses are known)
ALIGN_SKIP_BAND = 6.0      # cost of leaving a detected band unassigned (dust, well, dye front)
ALIGN_SKIP_SIZE = 5.0      # cost of a ladder size without band (missed, merged, run off)
ALIGN_CAP = 9.0            # one bad step never costs more than this
ALIGN_MAX_GAP = 3          # at most 3 bands or sizes skipped in a row
# assign_marker also runs the classic search and keeps it where its curve is smoother: by 10 %,
# or up to 1.5x rougher when it uses more bands (clean, evenly curved ladders, e.g. a merged
# pair right between its two sizes, or bands only a few pixels apart)
HYBRID_RATIO = 0.9
HYBRID_MORE = 1.5


def _align(ys: list[float], sizes: list[float], masses=None, vols=None) -> list[tuple[int, int]]:
    """Best monotone pairing (band index, size index) of sorted bands with the ladder.

    Dynamic programming over the last two pairs: each step costs the change of spacing per
    log(size) against the previous step and, with masses, the change of brightness per ng;
    unpaired bands and sizes cost a fixed penalty. About 0.1 s for 20 bands x 20 sizes."""
    m, n = len(ys), len(sizes)
    lg = [math.log10(float(s)) for s in sizes]
    lb = None
    if vols is not None and masses is not None and len(masses) == n:
        lb = [[math.log(max(float(v), 1e-9)) - math.log(max(float(w), 1e-9)) for w in masses] for v in vols]
    cap, g_ = ALIGN_CAP, ALIGN_MAX_GAP + 1
    inv_r, inv_b, wb = 1.0 / ALIGN_RATE_SD, 1.0 / ALIGN_BRIGHT_SD, ALIGN_BRIGHT_W
    sb, ss = ALIGN_SKIP_BAND, ALIGN_SKIP_SIZE

    def bright(i0, j0, i1, j1):
        if lb is None:
            return 0.0
        x = (lb[i1][j1] - lb[i0][j0]) * inv_b
        return wb * min(x * x, cap)

    # best[(i1, j1)][(i0, j0)] = (cost, pair before (i0, j0))
    best: dict = {}
    for i0 in range(m):
        for j0 in range(n):
            c0 = i0 * sb + j0 * ss
            for i1 in range(i0 + 1, min(m, i0 + g_ + 1)):
                for j1 in range(j0 + 1, min(n, j0 + g_ + 1)):
                    c = c0 + (i1 - i0 - 1) * sb + (j1 - j0 - 1) * ss + bright(i0, j0, i1, j1)
                    d = best.setdefault((i1, j1), {})
                    if (i0, j0) not in d or c < d[(i0, j0)][0]:
                        d[(i0, j0)] = (c, None)
    for i1 in range(m):
        for j1 in range(n):
            d = best.get((i1, j1))
            if not d:
                continue
            for (i0, j0), (c, _) in list(d.items()):
                r1 = (ys[i1] - ys[i0]) / (lg[j0] - lg[j1])
                if r1 <= 0:
                    continue
                for i2 in range(i1 + 1, min(m, i1 + g_ + 1)):
                    dy = ys[i2] - ys[i1]
                    if dy <= 0:
                        continue
                    for j2 in range(j1 + 1, min(n, j1 + g_ + 1)):
                        x = math.log(dy / (lg[j1] - lg[j2]) / r1) * inv_r
                        c2 = (c + min(x * x, cap) + (i2 - i1 - 1) * sb + (j2 - j1 - 1) * ss
                              + bright(i1, j1, i2, j2))
                        e = best.setdefault((i2, j2), {})
                        if (i1, j1) not in e or c2 < e[(i1, j1)][0]:
                            e[(i1, j1)] = (c2, (i0, j0))
    fin = None
    for (i1, j1), d in best.items():
        tail = (m - 1 - i1) * sb + (n - 1 - j1) * ss
        for prev, (c, _) in d.items():
            if fin is None or c + tail < fin[0]:
                fin = (c + tail, (i1, j1), prev)
    if fin is None:
        return []
    pairs, cur, prev = [fin[1]], fin[1], fin[2]
    while prev is not None:
        pairs.append(prev)
        cur, prev = prev, best[cur][prev][1]
    return pairs[::-1]


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
    if (vols is not None and masses is not None and len(vols) == len(masses) >= 4
            and np.std(vols) > 0 and np.std(masses) > 0):
        score += 0.05 * (1 - float(np.corrcoef(np.asarray(vols, float), np.asarray(masses, float))[0, 1]))
    return score


def assign_marker(band_ys: list[float], sizes: list[float], masses: list[float] | None = None,
                  volumes: list[float] | None = None, max_skip: int = 3, method: str = "align"
                  ) -> tuple[list[float], list[float], list[float] | None, list[str]]:
    """Pair detected marker bands with ladder sizes.

    ``align`` (default, 5+ bands and sizes): :func:`_align` pairs bands and sizes so that the
    spacing per log(size) - and, with masses, the brightness per ng - changes smoothly from step
    to step. Extra bands (dust, well, dye front) and missing sizes (missed, merged, run off) may
    occur together and anywhere.
    ``combos`` (and fewer than 5): equal counts -> direct pairing (top -> bottom); if counts
    differ by up to ``max_skip``, every way of dropping the surplus bands or sizes is tried and
    the smoothest size-vs-migration curve wins.
    A warning always says what was done so the user can check it.
    """
    order = np.argsort(band_ys)
    ys = [float(band_ys[i]) for i in order]
    vols = [float(volumes[i]) for i in order] if volumes is not None else None
    sizes = list(sizes)
    warnings: list[str] = []
    m, n = len(ys), len(sizes)
    if method == "align" and min(m, n) >= 5:
        pairs = _align(ys, sizes, masses, vols)
        if abs(m - n) <= max_skip:
            y_c, s_c, _, _ = assign_marker(ys, sizes, masses, vols, max_skip=max_skip, method="combos")
            iy, js = {y: i for i, y in enumerate(ys)}, {s: j for j, s in enumerate(sizes)}
            classic = [(iy[y], js[s]) for y, s in zip(y_c, s_c)]

            def rough(pp):
                return _fit_quality([ys[i] for i, _ in pp], [sizes[j] for _, j in pp],
                                    [vols[i] for i, _ in pp] if vols else None,
                                    [masses[j] for _, j in pp] if masses else None)
            if classic != pairs and len(classic) >= len(pairs):
                rc, ra = rough(classic), rough(pairs)
                if rc < HYBRID_RATIO * ra or (len(classic) > len(pairs) and rc < HYBRID_MORE * ra):
                    pairs = classic
        if len(pairs) >= 3:
            used_b = {i for i, _ in pairs}
            used_s = {j for _, j in pairs}
            dropped = [ys[i] for i in range(m) if i not in used_b]
            missing = [sizes[j] for j in range(n) if j not in used_s]
            if dropped:
                warnings.append(f"Marker: {m} bands detected, ladder has {n}. Ignored the band(s) at row(s) "
                                f"{', '.join(f'{y:.0f}' for y in dropped)} (best fit) - please check.")
            if missing:
                warnings.append(f"Marker: {len(pairs)} of {n} ladder bands found; assumed missing: "
                                f"{', '.join(f'{s:g}' for s in missing)} (best fit) - please check.")
            return ([ys[i] for i, _ in pairs], [sizes[j] for _, j in pairs],
                    [masses[j] for _, j in pairs] if masses else None, warnings)
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


def suggest_ladders(band_ys: list[float], volumes: list[float] | None = None, unit: str | None = None,
                    top: int = 3, max_skip: int = 2) -> list[dict]:
    """Which library ladders fit the bands of a marker lane best? (ranked, best first)

    Each ladder is paired with the bands like :func:`assign_marker` (up to ``max_skip`` extra or
    missing bands) and scored on (a) how smooth log(size) vs. migration is, (b) the band count and
    (c) the brightness pattern - datasheet masses, or ApE's relative intensities (e.g. the bright
    reference bands of many ladders). Ladders with identical sizes are reported once.
    Returns [{name, score, n_sizes, n_bands, matched, names_same_sizes}]."""
    ys = sorted(float(y) for y in band_ys)
    vols = None
    if volumes is not None:
        vols = [float(v) for _, v in sorted(zip(band_ys, volumes))]
    m = len(ys)
    if m < 4:
        return []
    seen: dict[tuple, dict] = {}
    for name, e in PRESETS.items():
        if unit and e.get("unit") != unit:
            continue
        sizes = list(e["sizes"])
        n = len(sizes)
        if abs(m - n) > max_skip or n < 4:
            continue
        ref = e.get("masses") or e.get("relative_intensity")
        y_a, s_a, r_a, _ = assign_marker(ys, sizes, ref, vols, max_skip=max_skip, method="combos")
        if len(s_a) < 4:
            continue
        shape = _fit_quality(y_a, s_a)
        score = shape + 0.03 * abs(m - n)
        if vols is not None:
            # brightness pattern: a ladder without one (no masses, all bands alike) gets a neutral
            # r = 0.5, otherwise ladders without data would beat ones whose pattern fits well
            r = 0.5
            v_a = [vols[ys.index(y)] for y in y_a]
            if r_a and np.std(r_a) > 0 and np.std(v_a) > 0:
                r = float(np.corrcoef(v_a, r_a)[0, 1])
            score += 0.1 * (1 - r)
        key = (e.get("unit"), tuple(sizes))
        hit = {"name": name, "score": round(float(score), 4), "n_sizes": n, "n_bands": m,
               "matched": len(s_a), "names_same_sizes": [name]}
        if key in seen:
            seen[key]["names_same_sizes"].append(name)
            # prefer the curated entry (it has masses) as the representative
            if PRESETS[seen[key]["name"]].get("library") == "ape" and e.get("library") != "ape":
                seen[key].update(name=name, score=hit["score"])
            continue
        seen[key] = hit
    return sorted(seen.values(), key=lambda h: h["score"])[:top]

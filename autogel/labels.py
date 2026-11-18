"""Band size labels: formatting, placement beside the lane, and zoom-aware decluttering.

Shared by the napari viewer and the exported annotated image so both look the same.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .markers import PRESETS


def fmt_size(size: float, unit: str = "bp") -> str:
    """10002 -> '10 kb', 1517 -> '1.5 kb', 900 -> '900', 37 kDa -> '37'."""
    if not np.isfinite(size):
        return ""
    if unit == "bp" and size >= 1000:
        kb = size / 1000
        return f"{kb:.0f} kb" if kb >= 9.95 else f"{kb:.1f} kb"
    return f"{size:.0f}"


@dataclass
class Label:
    band: int          # index into ga.bands
    lane: int
    y: float           # band peak row (data px)
    x: float           # anchor column (data px), just outside the lane
    side: str          # "left": text ends at x; "right": text starts at x
    text: str
    marker: bool
    priority: float    # higher = kept first when labels collide


def plan_labels(ga, df, lanes: set[int]) -> list[Label]:
    """One label per band of the given lanes, placed outside the lane.

    Marker lanes in the left half of the gel are labelled on their left, all other lanes on
    their right, so the labels of edge markers point away from the samples.
    """
    if df is None or not len(df) or not ga.lanes:
        return []
    sc = [c for c in df.columns if c.startswith("size_")]
    if not sc:
        return []
    unit = sc[0].split("_", 1)[1]
    sizes = {(int(r["lane"]) - 1, float(r["y_px"])): float(r[sc[0]]) for _, r in df.iterrows()}
    vols = {(int(r["lane"]) - 1, float(r["y_px"])): float(r["volume"]) for _, r in df.iterrows()}
    marker = {L - 1 for L in ga.params.marker_lanes}
    ref = set(PRESETS.get(ga.params.marker_preset or "", {}).get("reference_bands", []))
    mid = 0.5 * (ga.lanes[0].x + ga.lanes[-1].x)
    pad = max(3.0, 0.08 * float(np.median([l.width for l in ga.lanes])))
    out = []
    for k, b in enumerate(ga.bands):
        if b.lane not in lanes:
            continue
        key = (b.lane, round(b.y, 1))
        s = sizes.get(key)
        if s is None or not np.isfinite(s):
            continue
        l = ga.lanes[b.lane]
        is_marker = b.lane in marker
        left = is_marker and l.x < mid      # markers point outward, sample lanes label to the right
        prio = vols.get(key, 0.0)
        if is_marker and any(abs(s - r) / r < 0.02 for r in ref):
            prio += 1e18        # reference bands (e.g. 3 kb, 1 kb, 0.5 kb) always win
        o = l.dx(b.y)                       # curved lanes: label beside the lane at this row
        out.append(Label(band=k, lane=b.lane, y=b.y, x=(l.x0 + o - pad) if left else (l.x1 + o + pad),
                         side="left" if left else "right", text=fmt_size(s, unit),
                         marker=is_marker, priority=prio))
    return out


def declutter(labels: list[Label], min_gap: float) -> list[bool]:
    """Keep the highest-priority labels; drop any label closer than ``min_gap`` (data px)
    to an already kept label on the same lane side."""
    keep = [False] * len(labels)
    kept: dict[tuple[int, str], list[float]] = {}
    for i in sorted(range(len(labels)), key=lambda i: -labels[i].priority):
        lab = labels[i]
        ys = kept.setdefault((lab.lane, lab.side), [])
        if all(abs(lab.y - y) >= min_gap for y in ys):
            ys.append(lab.y)
            keep[i] = True
    return keep

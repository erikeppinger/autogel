"""Human-readable log messages and lab-book reports.

* ``describe(entry)``   one sentence per log entry (GUI status line, report)
* ``write_jsonl``       machine-readable log, one JSON object per line
* ``write_report_md``   Markdown report for electronic lab notebooks (eLabFTW, LabArchives,
                        OneNote, Benchling, Obsidian ... all accept Markdown or pasted Markdown)
"""
from __future__ import annotations

import hashlib
import json
import os
import platform

import numpy as np

from . import __version__


def _n(v, nd=1):
    return f"{v:.{nd}f}" if isinstance(v, (int, float)) else str(v)


def describe(e: dict) -> str:
    """Plain-language sentence for one log entry."""
    s, p, r = e.get("step", ""), e.get("params") or {}, e.get("result") or ""
    if s == "load":
        o = p.get("orientation", {})
        bits = []
        if o.get("rot90_k"):
            bits.append(f"rotated {90 * o['rot90_k']}°")
        if o.get("flip_lr") or o.get("flip_ud"):
            bits.append("flipped")
        if o.get("angle"):
            bits.append(f"deskewed {o['angle']:+.2f}°")
        if o.get("invert"):
            bits.append("inverted (dark bands)")
        shape = p.get("shape", ["?", "?"])
        msg = (f"Loaded {os.path.basename(p.get('path', ''))} ({shape[1]}×{shape[0]} px, "
               f"{p.get('dtype', '')}); orientation: {', '.join(bits) or 'unchanged (lanes already vertical)'}")
        return msg + (f". ⚠ {r}" if r and r != "ok" else "")
    if s == "orientation_change":
        return "Orientation changed by user"
    if s == "detect_lanes":
        return f"Detected {r}"
    if s == "confirm_lanes":
        return f"Lane count {p.get('n_lanes')}: {r}"
    if s == "set_lanes":
        return f"Lanes edited ({p.get('reason', '')}): {r}"
    if s == "set_lane_names":
        return f"Lane names set: {', '.join(p.get('names', []))}"
    if s == "detect_bands":
        lanes = p.get("lanes")
        where = "all lanes" if not lanes else f"lane(s) {', '.join(map(str, lanes))}"
        return (f"Band detection in {where}: {r} [background {p.get('bg_method')}, "
                f"radius {_n(p.get('bg_radius'))} px, min S/N {p.get('min_snr')}]")
    if s == "set_bands":
        return f"Bands edited ({p.get('reason', '')}): {r}"
    if s == "set_crop":
        return f"Analysed region ({p.get('reason', '')}): {r}"
    if s == "set_channel":
        return f"Analysed image: {r}"
    if s == "add_image":
        return (f"Added image {os.path.basename(p.get('path', ''))} as {r}"
                + (" (dark bands on light, inverted)" if p.get("invert") else ""))
    if s == "set_marker_channel":
        return f"Marker lanes read from: {r}"
    if s == "suggest_ladder":
        return (f"Ladder suggestion for lane {p.get('lane')}: {', '.join(p.get('ranking', []))}; {r}")
    if s == "set_marker":
        load = f", ladder loaded {p['load']} µg" if p.get("load") else ""
        return (f"Marker: {p.get('preset')} in lane(s) {', '.join(map(str, p.get('lanes', [])))}"
                f" ({len(p.get('sizes') or [])} bands{load})")
    if s == "set_columns":
        return f"Gels side by side: {r}" + (f" ({p['reason']})" if p.get("reason") else "")
    if s == "set_tiers":
        return f"Gel rows: {r}" + (f" ({p['reason']})" if p.get("reason") else "")
    if s == "set_size_range":
        return f"Valid size range: {r}" + (f" ({p['reason']})" if p.get("reason") else "")
    if s == "calibrate":
        return f"Size/mass calibration: {r}"
    if s == "quantify":
        return f"Quantified: {r}"
    if s == "calibrate+quantify":
        return f"Calibrated and quantified: {r}"
    if s == "batch":
        return f"Batch of {p.get('n_images')} gel(s) {r} -> {p.get('out_dir', '')}"
    if s == "export":
        return f"Exported to {p.get('out_dir', '')}"
    if s == "load_session":
        return f"Session restored: {r}"
    return f"{s}: {r}"


RECOMPUTE = ("calibrate", "quantify")


def condensed(log: list[dict]) -> list[tuple[dict, int]]:
    """Merge each run of automatic recomputations (calibrate/quantify after every edit) into a
    single entry: the last one of the run, with the number of recalculations."""
    out: list[tuple[dict, int]] = []
    run: list[dict] = []

    def flush():
        if not run:
            return
        n = max(1, sum(e["step"] == "quantify" for e in run))
        last_q = next((e for e in reversed(run) if e["step"] == "quantify"), run[-1])
        last_c = next((e for e in reversed(run) if e["step"] == "calibrate"), None)
        merged = dict(last_q)
        if last_c is not None and last_q["step"] == "quantify":
            merged = {**last_q, "step": "calibrate+quantify",
                      "result": f"{last_q['result']}; calibration: {last_c['result']}"}
        out.append((merged, n))
        run.clear()

    for e in log:
        if e["step"] in RECOMPUTE:
            run.append(e)
        else:
            flush()
            out.append((e, 1))
    flush()
    return out


def write_jsonl(ga, path: str):
    with open(path, "w", encoding="utf-8") as f:
        for e in ga.log:
            f.write(json.dumps({**e, "message": describe(e)}, ensure_ascii=False, default=str) + "\n")


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _md_table(rows: list[list], header: list[str]) -> str:
    def cell(v):
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            return ""
        return str(v).replace("|", "\\|")
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(cell(v) for v in r) + " |" for r in rows]
    return "\n".join(lines)


def write_report_md(ga, df, path: str, files: list[str] | None = None):
    """Lab-book report: provenance, settings, step log, per-lane results, warnings."""
    p = ga.params
    o = ga.gel.orientation
    try:
        digest = _sha256(ga.path)
    except OSError:
        digest = "n/a"
    started = ga.log[0]["time"] if ga.log else ""
    L = []
    L.append(f"# Gel analysis – {os.path.basename(ga.path)}\n")
    L.append(f"*autogel {__version__} · Python {platform.python_version()} · analysis started {started}*\n")

    L.append("## Image\n")
    L.append(_md_table([
        ["File", ga.path],
        ["SHA-256", f"`{digest}`"],
        ["Size / type", f"{ga.data.shape[1]} × {ga.data.shape[0]} px, {ga.gel.dtype}"],
        ["Saturation level", f"{ga.gel.saturation_level:.0f} ({int(ga.gel.raw_max.sum())} px at/above)"],
        ["Orientation", f"rot90×{o.get('rot90_k', 0)}, flip ↔ {bool(o.get('flip_lr'))}, "
                        f"flip ↕ {bool(o.get('flip_ud'))}, skew {o.get('angle', 0):+.2f}°, "
                        f"inverted {bool(o.get('invert'))}"],
    ], ["", ""]))

    L.append("\n## Settings\n")
    marker = ", ".join(map(str, p.marker_lanes)) or "–"
    L.append(_md_table([
        ["Lanes", f"{len(ga.lanes)} (requested: {p.n_lanes or 'auto'} per gel row), "
                  f"gel rows: {p.n_tiers or 1}, gels side by side: {p.n_columns or 1}"
                  + (f" split at row(s) {', '.join(f'{r:.0f}' for r in p.tier_rows)}"
                     if p.tier_rows and (p.n_tiers or 1) > 1 else "")],
        ["Background", f"{p.bg_method}, radius {p.bg_radius:.0f} px"],
        ["Band detection", f"min S/N {p.min_snr}, min rel. height {p.min_rel_height}, "
                           f"min height vs. gel {p.min_gel_rel}, "
                           f"min lane coverage {p.min_coverage}, smoothing σ {p.smooth_sigma} px"],
        ["Integration width", f"{p.lane_width_frac:.0%} of lane width"],
        ["Marker", f"{p.marker_preset or '–'} in lane(s) {marker}"],
        ["Valid size range", (f"{p.size_min or 0:g} – {p.size_max or '∞'} {p.marker_unit}"
                              if (p.size_min or p.size_max) else "automatic (whole ladder)")],
        ["Ladder load", f"{p.marker_load or 'as datasheet'} µg (datasheet: "
                        f"{p.marker_ref_load or 'library default'})"],
    ], ["Parameter", "Value"]))
    if ga.mass_fit:
        mf = ga.mass_fit
        how = ("signal per ng as a function of fragment size, fitted on the ladder"
               if mf.get("model") == "size" else f"volume = {mf['slope']:.4g} × {mf['unit']} (through origin)")
        L.append(f"\nMass curve: {how}; typical error of a ladder band left out of the fit "
                 f"±{mf.get('typical_error', float('nan')):.0%}, R² = {mf['r2']:.3f}, n = {mf['n']} "
                 f"unsaturated ladder bands, range {mf['min']:.0f}–{mf['max']:.0f} {mf['unit']}.\n")

    L.append("\n## Steps\n")
    rows = []
    for e, n in condensed(ga.log):
        rows.append([e["time"].replace("T", " "),
                     e["step"] + (f" (recalculated {n}×)" if n > 1 else ""), describe(e)])
    L.append(_md_table(rows, ["Time", "Step", "What happened"]))
    n_user = sum(1 for b in ga.bands if b.source != "auto")
    if n_user:
        L.append(f"\n{n_user} band(s) were added or edited by hand (`source = user` in the table).\n")

    L.append("\n## Results per lane\n")
    if df is not None and len(df):
        sc = [c for c in df.columns if c.startswith("size_")]
        mc = [c for c in df.columns if c.startswith("est_mass_")]
        rows = []
        for i, lane in enumerate(ga.lanes):
            sub = df[df.lane == i + 1]
            sizes = ", ".join(f"{v:.0f}" for v in sub[sc[0]] if np.isfinite(v)) if sc else ""
            sat = int(sub.saturated.sum()) if len(sub) else 0
            rows.append([i + 1, lane.name, "marker" if (i + 1) in p.marker_lanes else "",
                         len(sub), sizes,
                         f"{sub.volume.sum():.3g}" if len(sub) else "",
                         f"{sub[mc[0]].sum():.1f}" if mc and len(sub) else "",
                         sat or ""])
        unit = sc[0].split("_", 1)[1] if sc else ""
        L.append(_md_table(rows, ["Lane", "Name", "", "Bands", f"Sizes ({unit})", "Σ intensity (a.u.)",
                                  "Σ est. mass" + (f" ({mc[0].split('_')[-1]})" if mc else ""),
                                  "Saturated bands"]))
        L.append("\nFull band table: `*_bands.csv` / `*_autogel.xlsx`.\n")
    else:
        L.append("No bands quantified.\n")

    warn = list(ga.gel.warnings) + [w for c in ga.calibrations for w in c.warnings]
    if df is not None and len(df) and df.saturated.any():
        warn.append(f"{int(df.saturated.sum())} band(s) contain saturated pixels – their volume is "
                    "underestimated.")
    if df is not None and len(df) and "size_in_range" in df:
        out_r = int((df.size_in_range == False).sum())  # noqa: E712 (pandas)
        if out_r:
            warn.append(f"{out_r} band size(s) lie outside the valid size range / ladder "
                        "(extrapolated, less reliable).")
    L.append("\n## Warnings\n")
    L.append("\n".join(f"- ⚠ {w}" for w in warn) if warn else "None.")

    if files:
        L.append("\n\n## Files\n")
        L.append("\n".join(f"- `{os.path.basename(f)}`" for f in files))
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

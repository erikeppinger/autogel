"""Batch workflow: analyse many gels with the settings of a template (a finished gel/session).

Per gel: automatic orientation, lanes and bands with the template's settings, full export and a
session file (open it in the GUI to correct that gel and re-export). Combined outputs:

    batch_summary.csv   one row per gel (lanes, bands, warnings, mass-curve R², status)
    batch_bands.csv     all bands of all gels (column ``gel``)
    batch_autogel.xlsx  both tables
    batch_report.md     overview for the lab book, gels needing attention first
"""
from __future__ import annotations

import dataclasses
import datetime as _dt
import glob
import os
import traceback

import pandas as pd

from .analysis import GelAnalysis, Params
from .export import export_all

IMAGE_EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp")


def collect_images(items: list[str]) -> list[str]:
    """Files and folders (non-recursive) -> sorted list of image files."""
    out = []
    for it in items:
        if os.path.isdir(it):
            out += [p for p in glob.glob(os.path.join(it, "*")) if p.lower().endswith(IMAGE_EXT)]
        elif any(ch in it for ch in "*?"):
            out += [p for p in glob.glob(it) if p.lower().endswith(IMAGE_EXT)]
        elif os.path.isfile(it):
            out.append(it)
    seen, uniq = set(), []
    for p in sorted(out):
        k = os.path.normcase(os.path.abspath(p))
        if k not in seen:
            seen.add(k)
            uniq.append(p)
    return uniq


def template_params(source: Params | str) -> Params:
    """Settings to reuse: from Params or a session file. Gel-specific geometry is dropped."""
    if isinstance(source, str):
        import json
        with open(source, encoding="utf-8") as f:
            source = Params(**json.load(f)["params"])
    p = dataclasses.replace(source)
    p.tier_rows = None            # split rows are searched per image
    p.invert = None               # stain polarity detected per image
    p.sample_loads_ul = {}        # sample volumes belong to one gel
    return p


def run_batch(images: list[str], template: Params | str, out_dir: str, progress=None) -> pd.DataFrame:
    """Analyse ``images`` with the template's settings. ``progress(i, n, name)`` is called
    before each gel and may return False to cancel. Returns the summary table."""
    os.makedirs(out_dir, exist_ok=True)
    tp = template_params(template)
    rows, all_bands = [], []
    cancelled = False
    for i, path in enumerate(images):
        name = os.path.basename(path)
        if progress is not None and progress(i, len(images), name) is False:
            cancelled = True
            break
        row = {"gel": name, "path": os.path.abspath(path)}
        try:
            ga = GelAnalysis(path, dataclasses.replace(tp))
            ga.detect_lanes()
            ga.detect_bands()
            stem = os.path.splitext(name)[0]
            out = export_all(ga, os.path.join(out_dir, f"{stem}_autogel"))
            df = ga.quantify()
            warn = list(ga.gel.warnings) + [w for c in ga.calibrations for w in c.warnings]
            n_sat = int(df.saturated.sum()) if len(df) else 0
            expected = (tp.n_lanes or 0) * (tp.n_tiers or 1)
            if expected and len(ga.lanes) != expected:
                warn.append(f"{len(ga.lanes)} lanes found, {expected} expected")
            row.update({"status": "check" if warn else "ok", "lanes": len(ga.lanes), "bands": len(df),
                        "saturated_bands": n_sat,
                        "mass_r2": round(ga.mass_fit["r2"], 3) if ga.mass_fit else None,
                        "mass_typical_error": round(ga.mass_fit["typical_error"], 3) if ga.mass_fit else None,
                        "warnings": " | ".join(warn), "output": out,
                        "report": os.path.join(out, f"{stem}_report.md"),
                        "session": os.path.join(out, f"{stem}_session.json")})
            if len(df):
                all_bands.append(df.assign(gel=name)[["gel", *df.columns]])
        except Exception as e:  # one bad image must not stop the batch
            row.update({"status": "error", "warnings": f"{type(e).__name__}: {e}",
                        "traceback": traceback.format_exc(limit=3)})
        rows.append(row)
    summary = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["gel", "path", "status"])
    # how many of the chosen gels were analysed (the user may cancel in between)
    summary.attrs.update({"n_images": len(images), "cancelled": cancelled})
    bands = pd.concat(all_bands, ignore_index=True) if all_bands else pd.DataFrame()
    summary.to_csv(os.path.join(out_dir, "batch_summary.csv"), index=False)
    bands.to_csv(os.path.join(out_dir, "batch_bands.csv"), index=False)
    with pd.ExcelWriter(os.path.join(out_dir, "batch_autogel.xlsx")) as xw:
        summary.drop(columns=["traceback"], errors="ignore").to_excel(xw, sheet_name="summary", index=False)
        if len(bands):
            bands.to_excel(xw, sheet_name="bands", index=False)
    _write_batch_report(summary, tp, os.path.join(out_dir, "batch_report.md"))
    return summary


def _write_batch_report(summary: pd.DataFrame, tp: Params, path: str):
    from . import __version__
    order = {"error": 0, "check": 1, "ok": 2}
    s = summary.assign(_o=summary["status"].map(order)).sort_values(["_o", "gel"])
    lines = [f"# Batch analysis – {len(summary)} gel(s)\n",
             f"*autogel {__version__} · {_dt.datetime.now().isoformat(timespec='seconds')}*\n"]
    if summary.attrs.get("cancelled"):
        lines.append(f"**Cancelled** after {len(summary)} of {summary.attrs.get('n_images')} gel(s); "
                     "the remaining gels were not analysed.\n")
    lines += ["## Settings used for all gels\n",
             f"- lanes per gel row: {tp.n_lanes or 'auto'}; gel rows: {tp.n_tiers or 1}",
             f"- background: {tp.bg_method}; min S/N {tp.min_snr}; min lane coverage {tp.min_coverage}",
             f"- marker: {tp.marker_preset or '–'} in lane(s) {', '.join(map(str, tp.marker_lanes)) or '–'}",
             f"- valid size range: {tp.size_min or 0:g} – {tp.size_max or '∞'}"
             if (tp.size_min or tp.size_max) else "- valid size range: automatic",
             "", "## Gels (needing attention first)\n",
             "| Gel | Status | Lanes | Bands | Saturated | Mass R² | Notes |", "|---|---|---|---|---|---|---|"]
    icon = {"ok": "✔ ok", "check": "⚠ check", "error": "✖ error"}

    def cell(row, k):
        x = row.get(k)
        if x is None or (isinstance(x, float) and pd.isna(x)):
            return ""
        if isinstance(x, float) and x.is_integer() and k != "mass_r2":
            return str(int(x))
        return str(x)

    for _, r in s.iterrows():
        note = cell(r, "warnings").replace("|", "/")
        lines.append(f"| {r['gel']} | {icon.get(r['status'], r['status'])} | {cell(r, 'lanes')} "
                     f"| {cell(r, 'bands')} | {cell(r, 'saturated_bands')} | {cell(r, 'mass_r2')} | {note} |")
    lines += ["", "Each gel has its own folder `<gel>_autogel/` with report, tables and a session file. "
              "To correct a gel, open its session in the GUI, fix lanes/bands and export again."]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

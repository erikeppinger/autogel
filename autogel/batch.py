"""Batch workflow: analyse many gels with the settings of a template (a finished gel/session).

Per gel: automatic orientation, lanes and bands with the template's settings, full export and a
session file (open it in the GUI to correct that gel and re-export). Combined outputs:

    batch_summary.csv   one row per gel (lanes, bands, warnings, mass-curve R², status)
    batch_bands.csv     all bands of all gels (column ``gel``)
    batch_compare.csv   with a reference lane: every sample's bands relative to the reference
                        sample, matched by size across gels (same reference sample on every gel)
    batch_autogel.xlsx  all tables
    batch_report.md     overview for the lab book, gels needing attention first
"""
from __future__ import annotations

import dataclasses
import datetime as _dt
import glob
import os
import traceback

import numpy as np
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


# ------------------------------------------------------------------ sample sheet
# One row per gel; empty cells keep the template's setting. Header spellings are forgiving.
SHEET_COLUMNS = ["image", "lanes", "gel_rows", "marker_lanes", "ladder", "lane_names", "reference_lane"]
_ALIASES = {"image": "image", "picture": "image", "file": "image", "gel": "image",
            "lanes": "lanes", "number_of_lanes": "lanes", "n_lanes": "lanes",
            "gel_rows": "gel_rows", "rows": "gel_rows",
            "marker_lanes": "marker_lanes", "marker": "marker_lanes", "ladder_lanes": "marker_lanes",
            "ladder": "ladder", "marker_name": "ladder",
            "lane_names": "lane_names", "names": "lane_names", "samples": "lane_names", "sample_names": "lane_names",
            "reference_lane": "reference_lane", "reference": "reference_lane", "control_lane": "reference_lane"}


def _cell(x) -> str:
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return ""
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x).strip()


def _split(text: str) -> list[str]:
    """'1, 12' / '1; 12' -> ['1', '12'] (a cell may use either separator)."""
    return [t.strip() for t in text.replace(";", ",").split(",") if t.strip()]


def read_sample_sheet(path: str) -> list[dict]:
    """Sample sheet (.xlsx/.xls/.csv) -> one dict per gel:
    {"row", "path", "n_lanes", "n_tiers", "marker_lanes", "ladder", "lane_names", "reference_lane"}
    (None = use the template's setting).
    Picture paths may be relative to the sheet's folder. Raises ValueError with a plain message."""
    from .markers import PRESETS
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, dtype=object)
    else:   # Excel writes ';' in many countries, ',' in others: let pandas detect it
        df = pd.read_csv(path, sep=None, engine="python", dtype=object, encoding="utf-8-sig")
    cols = {}
    for c in df.columns:
        key = _ALIASES.get(str(c).strip().lower().replace(" ", "_").replace("-", "_"))
        if key and key not in cols.values():
            cols[c] = key
    if "image" not in cols.values():
        raise ValueError(f"The sample sheet needs a column 'image' with the picture file names "
                         f"(found: {', '.join(map(str, df.columns))}).")
    df = df[list(cols)].rename(columns=cols)
    base = os.path.dirname(os.path.abspath(path))
    out, problems = [], []
    for k, r in df.iterrows():
        row = k + 2                              # as the user sees it in Excel (header = row 1)
        img = _cell(r.get("image"))
        if not img:
            continue
        g = {"row": row, "path": img if os.path.isabs(img) else os.path.join(base, img),
             "n_lanes": None, "n_tiers": None, "marker_lanes": None, "ladder": None, "lane_names": None,
             "reference_lane": None}
        try:
            if _cell(r.get("lanes")):
                g["n_lanes"] = int(float(_cell(r.get("lanes"))))
            if _cell(r.get("gel_rows")):
                g["n_tiers"] = int(float(_cell(r.get("gel_rows"))))
            if _cell(r.get("reference_lane")):
                g["reference_lane"] = int(float(_cell(r.get("reference_lane"))))
            if _cell(r.get("marker_lanes")):
                g["marker_lanes"] = [int(float(t)) for t in _split(_cell(r.get("marker_lanes")))]
        except ValueError:
            problems.append(f"row {row}: lanes, gel rows, marker lanes and reference lane must be numbers")
            continue
        lad = _cell(r.get("ladder"))
        if lad:
            match = [n for n in PRESETS if n.lower() == lad.lower()]
            if not match:
                problems.append(f"row {row}: ladder '{lad}' is not in autogel's list - copy its name "
                                "exactly from the ladder list in the panel")
                continue
            g["ladder"] = match[0]
        if _cell(r.get("lane_names")):
            g["lane_names"] = _split(_cell(r.get("lane_names")))
        if g["n_lanes"] is not None and g["n_lanes"] < 1:
            problems.append(f"row {row}: number of lanes must be at least 1")
            continue
        if g["marker_lanes"] and g["n_lanes"] and max(g["marker_lanes"]) > g["n_lanes"] * (g["n_tiers"] or 1):
            problems.append(f"row {row}: marker lane {max(g['marker_lanes'])} is beyond the {g['n_lanes']} lanes")
            continue
        if g["reference_lane"] and g["n_lanes"] and g["reference_lane"] > g["n_lanes"] * (g["n_tiers"] or 1):
            problems.append(f"row {row}: reference lane {g['reference_lane']} is beyond the {g['n_lanes']} lanes")
            continue
        out.append(g)
    if problems:
        raise ValueError("Please correct the sample sheet:\n" + "\n".join(problems))
    if not out:
        raise ValueError("The sample sheet lists no pictures (column 'image' is empty).")
    return out


def write_sample_sheet(path: str, images: list[str], template: Params | None = None) -> str:
    """Excel sheet to fill in: one row per picture, pre-filled with the template's settings."""
    tp = template or Params()
    base = os.path.dirname(os.path.abspath(path))
    rows = []
    for p in images:
        rel = os.path.relpath(os.path.abspath(p), base)
        rows.append({"image": p if rel.startswith("..") else rel,
                     "lanes": tp.n_lanes or "", "gel_rows": tp.n_tiers or 1,
                     "marker_lanes": ", ".join(map(str, tp.marker_lanes)),
                     "ladder": tp.marker_preset if tp.marker_preset and tp.marker_preset != "custom" else "",
                     "lane_names": "", "reference_lane": tp.reference_lane or ""})
    df = pd.DataFrame(rows, columns=SHEET_COLUMNS)
    with pd.ExcelWriter(path) as xw:
        df.to_excel(xw, sheet_name="gels", index=False)
        help_rows = [("image", "picture file (next to this sheet, or a full path)"),
                     ("lanes", "number of lanes per gel row, ladder lanes included"),
                     ("gel_rows", "gels stacked in one picture (usually 1)"),
                     ("marker_lanes", "lane(s) with the ladder, e.g. 1, 12"),
                     ("ladder", "ladder name exactly as in autogel's ladder list"),
                     ("lane_names", "names from left to right, separated by commas"),
                     ("reference_lane", "lane with the reference sample: every band is shown relative to "
                                        "the same band there - load the same sample on every gel to compare gels"),
                     ("", "Empty cells use the settings of the gel open in autogel.")]
        pd.DataFrame(help_rows, columns=["column", "what to enter"]).to_excel(xw, sheet_name="how to fill in",
                                                                             index=False)
    return path


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


def _gel_params(tp: Params, g: dict | None) -> Params:
    """The template's settings with this gel's sample-sheet entries applied."""
    p = dataclasses.replace(tp)
    if not g:
        return p
    if g.get("n_lanes") is not None:
        p.n_lanes = g["n_lanes"]
    if g.get("n_tiers") is not None:
        p.n_tiers = g["n_tiers"]
    if g.get("marker_lanes") is not None:
        p.marker_lanes = list(g["marker_lanes"])
    if g.get("reference_lane") is not None:
        p.reference_lane = g["reference_lane"]
    if g.get("ladder"):
        from .markers import PRESETS
        pre = PRESETS[g["ladder"]]
        p.marker_preset, p.marker_sizes, p.marker_unit = g["ladder"], list(pre["sizes"]), pre["unit"]
        p.marker_masses = list(pre["masses"]) if pre.get("masses") else None
    return p


def run_batch(images: list[str], template: Params | str, out_dir: str, progress=None,
              sheet: list[dict] | None = None) -> pd.DataFrame:
    """Analyse ``images`` with the template's settings. ``progress(i, n, name)`` is called
    before each gel and may return False to cancel. ``sheet`` (from ``read_sample_sheet``) gives
    gels their own lanes, marker lanes, ladder and lane names; its pictures are analysed after
    ``images``. Returns the summary table."""
    os.makedirs(out_dir, exist_ok=True)
    tp = template_params(template)
    jobs = [(p, None) for p in images] + [(g["path"], g) for g in (sheet or [])]
    rows, all_bands = [], []
    cancelled = False
    for i, (path, g) in enumerate(jobs):
        name = os.path.basename(path)
        if progress is not None and progress(i, len(jobs), name) is False:
            cancelled = True
            break
        row = {"gel": name, "path": os.path.abspath(path)}
        if g is not None:
            row["sheet_row"] = g["row"]
        try:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"picture not found: {path}")
            p = _gel_params(tp, g)
            ga = GelAnalysis(path, p)
            ga.detect_lanes()
            ga.detect_bands()
            warn = []
            if g is not None and g.get("lane_names"):
                names = g["lane_names"]
                if len(names) != len(ga.lanes):
                    warn.append(f"{len(names)} lane names for {len(ga.lanes)} lanes")
                ga.set_lane_names(names[: len(ga.lanes)])
            stem = os.path.splitext(name)[0]
            out = export_all(ga, os.path.join(out_dir, f"{stem}_autogel"))
            df = ga.quantify()
            warn += list(ga.gel.warnings) + [w for c in ga.calibrations for w in c.warnings]
            n_sat = int(df.saturated.sum()) if len(df) else 0
            expected = (p.n_lanes or 0) * (p.n_tiers or 1)
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
    summary.attrs.update({"n_images": len(jobs), "cancelled": cancelled, "sheet": bool(sheet)})
    bands = pd.concat(all_bands, ignore_index=True) if all_bands else pd.DataFrame()
    across = compare_across_gels(bands, tp.group_tol)
    summary.attrs["compare"] = bool(len(across))
    if len(across):
        across.to_csv(os.path.join(out_dir, "batch_compare.csv"), index=False)
    summary.to_csv(os.path.join(out_dir, "batch_summary.csv"), index=False)
    bands.to_csv(os.path.join(out_dir, "batch_bands.csv"), index=False)
    with pd.ExcelWriter(os.path.join(out_dir, "batch_autogel.xlsx")) as xw:
        summary.drop(columns=["traceback"], errors="ignore").to_excel(xw, sheet_name="summary", index=False)
        if len(bands):
            bands.to_excel(xw, sheet_name="bands", index=False)
        if len(across):
            across.to_excel(xw, sheet_name="compare across gels", index=False)
    _write_batch_report(summary, tp, os.path.join(out_dir, "batch_report.md"))
    return summary


def compare_across_gels(bands: pd.DataFrame, tol: float = 0.04) -> pd.DataFrame:
    """Samples (gel, lane) x band sizes table of ``fold_change`` - each band relative to the same band
    in its gel's reference lane. With the same reference sample on every gel, these values can be
    compared between gels. Bands are matched across gels by size (within ``tol``); needs a ladder."""
    need = {"fold_change", "group_size", "gel", "lane", "lane_name", "is_marker"}
    if not len(bands) or not need <= set(bands.columns):
        return pd.DataFrame()
    b = bands[~bands["is_marker"].astype(bool) & bands["fold_change"].notna() & bands["group_size"].notna()]
    if not len(b):
        return pd.DataFrame()
    b = b.assign(_k=-np.log10(b["group_size"].astype(float))).sort_values("_k")
    cluster, start, c = [], None, -1
    for k in b["_k"]:
        if start is None or k - start > np.log10(1 + tol):
            start, c = k, c + 1
        cluster.append(c)
    b = b.assign(_c=cluster)
    label = b.groupby("_c")["group_size"].median().map(lambda v: f"{float(f'{v:.3g}'):g}")   # 2500, 799
    b = b.assign(band=b["_c"].map(label))
    m = b.pivot_table(index=["gel", "lane", "lane_name"], columns="band", values="fold_change", aggfunc="first")
    m = m.reindex(columns=list(dict.fromkeys(b.sort_values("_k")["band"])))
    m.columns = [f"{c} (relative to reference)" for c in m.columns]
    return m.reset_index()


def _write_batch_report(summary: pd.DataFrame, tp: Params, path: str):
    from . import __version__
    order = {"error": 0, "check": 1, "ok": 2}
    s = summary.assign(_o=summary["status"].map(order)).sort_values(["_o", "gel"])
    lines = [f"# Batch analysis – {len(summary)} gel(s)\n",
             f"*autogel {__version__} · {_dt.datetime.now().isoformat(timespec='seconds')}*\n"]
    if summary.attrs.get("cancelled"):
        lines.append(f"**Cancelled** after {len(summary)} of {summary.attrs.get('n_images')} gel(s); "
                     "the remaining gels were not analysed.\n")
    if summary.attrs.get("sheet"):
        lines.append("Gels from the **sample sheet** use their own lanes, gel rows, marker lanes, ladder and "
                     "lane names where the sheet gives them; empty cells use the settings below.\n")
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
    if summary.attrs.get("compare"):
        lines += ["", "**Comparing gels:** `batch_compare.csv` shows every sample's bands relative to the same "
                  "band in its gel's reference lane. This is comparable between gels only if the same "
                  "reference sample was loaded on every gel."]
    lines += ["", "Each gel has its own folder `<gel>_autogel/` with report, tables and a session file. "
              "To correct a gel, open its session in the GUI, fix lanes/bands and export again."]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

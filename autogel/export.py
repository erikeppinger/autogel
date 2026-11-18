"""Export: tables (CSV/XLSX), lane profiles (CSV + PDF), annotated image, log, session."""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg", force=False)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from .analysis import GelAnalysis  # noqa: E402


def _display(a: np.ndarray) -> tuple[float, float]:
    return float(np.percentile(a, 1)), float(np.percentile(a, 99.7))


def profiles_long(ga: GelAnalysis) -> pd.DataFrame:
    rows = []
    for i, l in enumerate(ga.lanes):
        pr = ga.lane_profile(i)
        rows.append(pd.DataFrame({"lane": i + 1, "lane_name": l.name, "y_px": pr["y"].astype(int),
                                  "raw": pr["raw"], "baseline": pr["baseline"], "corrected": pr["corr"],
                                  "saturated_px": pr["sat"]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_lane(ga: GelAnalysis, i: int, ax, df: pd.DataFrame | None = None):
    pr = ga.lane_profile(i)
    y = pr["y"]
    ax.plot(y, pr["raw"], color="0.6", lw=0.8, label="raw")
    ax.plot(y, pr["baseline"], color="#E69F00", lw=0.8, ls="--", label="background")
    ax.plot(y, pr["corr"], color="k", lw=1.0, label="corrected")
    sat = pr["sat"] > 0
    if sat.any():
        ax.plot(y[sat], pr["corr"][sat], ".", color="#D55E00", ms=2, label="saturated")
    bands = [b for b in ga.bands if b.lane == i]
    sub = df[df["lane"] == i + 1] if df is not None and len(df) else None
    for k, b in enumerate(bands):
        sel = (y >= b.y0) & (y <= b.y1)
        ax.fill_between(y[sel], 0, pr["corr"][sel], color="#56B4E9", alpha=0.4)
        lab = f"{k + 1}"
        if sub is not None and k < len(sub):
            sc = [c for c in sub.columns if c.startswith("size_")]
            if sc and np.isfinite(sub.iloc[k][sc[0]]):
                lab += f"\n{sub.iloc[k][sc[0]]:.0f}"
        ax.annotate(lab, (b.y, pr["corr"][sel].max() if sel.any() else 0), ha="center",
                    va="bottom", fontsize=6)
    ax.set_xlabel("migration (px from top of oriented image)")
    ax.set_ylabel("integrated intensity")
    ax.set_title(f"Lane {i + 1}: {ga.lanes[i].name}", fontsize=9)
    ax.legend(fontsize=6, loc="upper right")


def annotated_image(ga: GelAnalysis, path: str, df: pd.DataFrame | None = None):
    """Full-resolution PNG with lanes, bands and sizes, drawn with Pillow on an 8-bit
    copy (matplotlib's image pipeline needed ~0.8 GB for a 5 MP gel; this needs ~30 MB)."""
    from PIL import Image, ImageDraw, ImageFont
    a = ga.data
    lo, hi = _display(a)
    g = np.clip((a - lo) * np.float32(255.0 / max(hi - lo, 1e-9)), 0, 255).astype(np.uint8)
    rgb = np.repeat(g[:, :, None], 3, axis=2)
    del g
    rgb[ga.gel.raw_max] = (0xD5, 0x5E, 0x00)          # Okabe-Ito vermillion = saturated
    im = Image.fromarray(rgb)
    del rgb
    d = ImageDraw.Draw(im)
    fsize = max(10, int(np.median([l.width for l in ga.lanes]) / 3)) if ga.lanes else 12
    try:
        font = ImageFont.load_default(size=fsize)
        small = ImageFont.load_default(size=max(8, int(fsize * 0.7)))
    except TypeError:                                   # Pillow < 10.1
        font = small = ImageFont.load_default()
    marker = {L - 1 for L in ga.params.marker_lanes}
    for i, l in enumerate(ga.lanes):
        d.rectangle([l.x0, min(l.y0, l.y1), l.x1, max(l.y0, l.y1)],
                    outline="#E69F00" if i in marker else "#56B4E9", width=2)
        d.text((l.x, max(0, l.y0 - 4)), l.name, fill="white", font=font, anchor="ms",
               stroke_width=2, stroke_fill="black")
    sizes = {}
    if df is not None and len(df):
        sc = [c for c in df.columns if c.startswith("size_")]
        if sc:
            sizes = {(r.lane - 1, r.y_px): r[sc[0]] for _, r in df.iterrows()}
    for b in ga.bands:
        l = ga.lanes[b.lane]
        d.rectangle([l.x0, b.y0, l.x1, b.y1 + 1],
                    outline="#F0E442" if b.source == "auto" else "#CC79A7", width=1)
        s = sizes.get((b.lane, round(b.y, 1)))
        if s is not None and np.isfinite(s):
            d.text((l.x1 + 2, b.y), f"{s:.0f}", fill="#F0E442", font=small, anchor="lm",
                   stroke_width=2, stroke_fill="black")
    im.save(path, optimize=False)


def write_log_txt(ga: GelAnalysis, path: str):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"autogel analysis log\nimage: {ga.path}\n\n")
        for e in ga.log:
            f.write(f"[{e['time']}] {e['step']}\n")
            if e["params"]:
                f.write(f"    params: {json.dumps(e['params'], default=str)}\n")
            if e["result"]:
                f.write(f"    result: {e['result']}\n")
        if ga.gel.warnings:
            f.write("\nWARNINGS\n" + "\n".join(" - " + w for w in ga.gel.warnings) + "\n")
        for c in ga.calibrations:
            for w in c.warnings:
                f.write(" - " + w + "\n")


def export_all(ga: GelAnalysis, out_dir: str | None = None) -> str:
    stem = os.path.splitext(os.path.basename(ga.path))[0]
    out_dir = out_dir or os.path.join(os.path.dirname(ga.path), f"{stem}_autogel")
    os.makedirs(out_dir, exist_ok=True)
    df = ga.quantify()
    lanes = ga.lanes_table()
    prof = profiles_long(ga)
    ga._log("export", {"out_dir": os.path.abspath(out_dir)})

    df.to_csv(os.path.join(out_dir, f"{stem}_bands.csv"), index=False)
    lanes.to_csv(os.path.join(out_dir, f"{stem}_lanes.csv"), index=False)
    prof.to_csv(os.path.join(out_dir, f"{stem}_profiles.csv"), index=False)
    cal_rows = []
    for c in ga.calibrations:
        for k, (y, s) in enumerate(zip(c.ys, c.sizes)):
            cal_rows.append({"marker_lane": c.lane + 1, "y_px": y, "size": s,
                             "mass_datasheet": c.masses[k] if c.masses else None,
                             "volume": c.volumes[k] if c.volumes else None})
    cal = pd.DataFrame(cal_rows)
    log = pd.DataFrame([{**e, "params": json.dumps(e["params"], default=str)} for e in ga.log])
    with pd.ExcelWriter(os.path.join(out_dir, f"{stem}_autogel.xlsx")) as xw:
        df.to_excel(xw, sheet_name="bands", index=False)
        lanes.to_excel(xw, sheet_name="lanes", index=False)
        cal.to_excel(xw, sheet_name="calibration", index=False)
        if ga.mass_fit:
            pd.DataFrame([{k: v for k, v in ga.mass_fit.items() if k != "points"}]).to_excel(
                xw, sheet_name="mass_fit", index=False)
        log.to_excel(xw, sheet_name="log", index=False)

    with PdfPages(os.path.join(out_dir, f"{stem}_profiles.pdf")) as pdf:
        for i in range(len(ga.lanes)):
            fig, ax = plt.subplots(figsize=(8, 3.5))
            plot_lane(ga, i, ax, df)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        if ga.calibrations:
            fig, axs = plt.subplots(1, 2 if ga.mass_fit else 1, figsize=(9, 3.5), squeeze=False)
            for c in ga.calibrations:
                axs[0, 0].semilogy(c.ys, c.sizes, "o-", ms=3, label=f"marker lane {c.lane + 1}")
            axs[0, 0].set_xlabel("y (px)")
            axs[0, 0].set_ylabel(f"size ({ga.params.marker_unit})")
            axs[0, 0].legend(fontsize=7)
            if ga.mass_fit:
                pts = np.array(ga.mass_fit["points"])
                axs[0, 1].plot(pts[:, 0], pts[:, 1], "o", ms=3)
                xx = np.linspace(0, pts[:, 0].max() * 1.05, 10)
                axs[0, 1].plot(xx, ga.mass_fit["slope"] * xx, "-")
                axs[0, 1].set_xlabel(f"ladder band mass ({ga.mass_fit['unit']})")
                axs[0, 1].set_ylabel("volume")
                axs[0, 1].set_title(f"R² = {ga.mass_fit['r2']:.3f}", fontsize=9)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    annotated_image(ga, os.path.join(out_dir, f"{stem}_annotated.png"), df)
    write_log_txt(ga, os.path.join(out_dir, f"{stem}_log.txt"))
    ga.save_session(os.path.join(out_dir, f"{stem}_session.json"))
    return out_dir

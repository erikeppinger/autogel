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
from .compare import group_matrix  # noqa: E402
from .labels import declutter, plan_labels  # noqa: E402
from .report import describe, write_jsonl, write_report_md  # noqa: E402


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
        if l.bend:                                      # curved lane: polygon along its centre line
            d.polygon([(float(x), float(y)) for y, x in l.outline()],
                      outline="#E69F00" if i in marker else "#56B4E9", width=2)
        else:
            d.rectangle([l.x0, min(l.y0, l.y1), l.x1, max(l.y0, l.y1)],
                        outline="#E69F00" if i in marker else "#56B4E9", width=2)
        d.text((l.x, max(0, l.y0 - 4)), l.name, fill="white", font=font, anchor="ms",
               stroke_width=2, stroke_fill="black")
    for b in ga.bands:
        l = ga.lanes[b.lane]
        o = l.dx(b.y)
        d.rectangle([l.x0 + o, b.y0, l.x1 + o, b.y1 + 1],
                    outline="#F0E442" if b.source == "auto" else "#CC79A7", width=1)
    # size labels beside the lanes (markers outward, samples to the right), decluttered
    labels = plan_labels(ga, df, set(range(len(ga.lanes))))
    line_h = small.getbbox("0")[3] + 3 if hasattr(small, "getbbox") else 14
    for lab, keep in zip(labels, declutter(labels, line_h)):
        if keep:
            d.text((lab.x, lab.y), lab.text, fill="#E69F00" if lab.marker else "white", font=small,
                   anchor="rm" if lab.side == "left" else "lm", stroke_width=2, stroke_fill="black")
    im.save(path, optimize=False)


def write_log_txt(ga: GelAnalysis, path: str):
    """Plain-text log: one readable line per step, raw parameters indented below."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"autogel analysis log\nimage: {ga.path}\n\n")
        for e in ga.log:
            f.write(f"[{e['time']}] {e['step']}: {describe(e)}\n")
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
    best = "fold_change" if "fold_change" in df else "normalised" if "normalised" in df else "volume"
    m = group_matrix(df, best)
    if len(m):
        m.to_csv(os.path.join(out_dir, f"{stem}_compare_{best}.csv"), index=False)
    lanes.to_csv(os.path.join(out_dir, f"{stem}_lanes.csv"), index=False)
    prof.to_csv(os.path.join(out_dir, f"{stem}_profiles.csv"), index=False)
    cal_rows = []
    for c in ga.calibrations:
        for k, (y, s) in enumerate(zip(c.ys, c.sizes)):
            cal_rows.append({"marker_lane": c.lane + 1, "y_px": y, "size": s,
                             "mass_datasheet": c.masses[k] if c.masses else None,
                             "volume": c.volumes[k] if c.volumes else None})
    cal = pd.DataFrame(cal_rows)
    log = pd.DataFrame([{"time": e["time"], "step": e["step"], "message": describe(e),
                         "result": e["result"], "params": json.dumps(e["params"], default=str)}
                        for e in ga.log])
    with pd.ExcelWriter(os.path.join(out_dir, f"{stem}_autogel.xlsx")) as xw:
        df.to_excel(xw, sheet_name="bands", index=False)
        lanes.to_excel(xw, sheet_name="lanes", index=False)
        cal.to_excel(xw, sheet_name="calibration", index=False)
        if ga.mass_fit:
            pd.DataFrame([{k: v for k, v in ga.mass_fit.items() if k != "points"}]).to_excel(
                xw, sheet_name="mass_fit", index=False)
        # band groups x lanes: the same band compared across lanes
        for value in ("volume", "normalised", "fold_change"):
            m = group_matrix(df, value)
            if len(m):
                m.to_excel(xw, sheet_name=f"compare_{value}"[:31], index=False)
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
                mf = ga.mass_fit
                pts = np.array(mf["points"])
                ax = axs[0, 1]
                if mf["model"] == "size":           # signal per ng vs fragment size
                    ax.semilogx(pts[:, 2], pts[:, 1] / pts[:, 0], "o", ms=3)
                    ll = np.linspace(*mf["log_size_range"], 100)
                    ax.semilogx(10 ** ll, 10 ** np.polyval(mf["coef"], ll), "-")
                    ax.set_xlabel(f"ladder band size ({ga.params.marker_unit})")
                    ax.set_ylabel(f"intensity per {mf['unit']}")
                else:
                    ax.plot(pts[:, 0], pts[:, 1], "o", ms=3)
                    xx = np.linspace(0, pts[:, 0].max() * 1.05, 10)
                    ax.plot(xx, mf["slope"] * xx, "-")
                    ax.set_xlabel(f"ladder band mass ({mf['unit']})")
                    ax.set_ylabel("integrated intensity (a.u.)")
                ax.set_title(f"mass curve: typical error ±{mf['typical_error']:.0%} (leave-one-out)", fontsize=9)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    annotated_image(ga, os.path.join(out_dir, f"{stem}_annotated.png"), df)
    write_log_txt(ga, os.path.join(out_dir, f"{stem}_log.txt"))
    write_jsonl(ga, os.path.join(out_dir, f"{stem}_log.jsonl"))
    ga.save_session(os.path.join(out_dir, f"{stem}_session.json"))
    files = sorted(f for f in os.listdir(out_dir) if f.startswith(stem + "_") and not f.endswith("_report.md"))
    write_report_md(ga, df, os.path.join(out_dir, f"{stem}_report.md"),
                    files=[f"{stem}_report.md", *files])
    return out_dir

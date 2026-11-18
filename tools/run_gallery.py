"""Render the simulated gel gallery, analyse every gel and write a visual summary.

    python tools/run_gallery.py OUTPUT_DIR

Writes OUTPUT_DIR/summary.md (scores per scenario), OUTPUT_DIR/overview.png (annotated
thumbnails) and one full autogel export per scenario.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

from tests.gallery import SCENARIOS, evaluate  # noqa: E402


def main(out: str):
    os.makedirs(out, exist_ok=True)
    rows = []
    for sc in SCENARIOS:
        r = evaluate(sc, out, export=True)
        rows.append(r)
        print(f"{sc.name:28s} lanes {r['lanes_found']}/{sc.n_lanes}  bands found {r['bands_found']} "
              f"missed {r['bands_missed']} extra {r['extra_bands']}  size err {r['median_size_error_pct']} %")
    with open(os.path.join(out, "summary.md"), "w", encoding="utf-8") as f:
        f.write("# autogel – simulated gel gallery\n\n")
        f.write("| Scenario | Conditions | Lanes | Ladder bands | Sample bands found / missed / extra "
                "| Size error median / max (%) | Ratio error (%) |\n|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['scenario']} | {r['description']} | {r['lanes_found']}/{r['lanes_expected']} "
                    f"| {r['ladder_bands']} | {r['bands_found']} / {r['bands_missed']} / {r['extra_bands']} "
                    f"| {r['median_size_error_pct']} / {r['max_size_error_pct']} "
                    f"| {r['median_ratio_error_pct']} |\n")
    n = len(rows)
    cols = 5
    fig, axs = plt.subplots((n + cols - 1) // cols, cols, figsize=(cols * 3.4, ((n + cols - 1) // cols) * 4))
    for ax in axs.ravel():
        ax.set_axis_off()
    for ax, r in zip(axs.ravel(), rows):
        png = [f for f in os.listdir(r["out"]) if f.endswith("_annotated.png")][0]
        im = Image.open(os.path.join(r["out"], png))
        im.thumbnail((700, 700))
        ax.imshow(im)
        ax.set_title(f"{r['scenario']}\nsize err {r['median_size_error_pct']} %, extra {r['extra_bands']}",
                     fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "overview.png"), dpi=110)
    print("written:", os.path.join(out, "summary.md"), os.path.join(out, "overview.png"))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "gallery_out")

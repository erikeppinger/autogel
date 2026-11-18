"""Comparing lanes: band groups, normalisation (Western blots), RNA integrity.

All functions work on the band table from :meth:`GelAnalysis.quantify` (plus the analysis for
profiles) and only ADD columns, so the plain densitometry results stay unchanged.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd


def channel_column(name: str) -> str:
    """Band-table column holding a channel's volume, e.g. 'volume_page_2'."""
    return "volume_" + re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()


def add_band_groups(df: pd.DataFrame, size_col: str | None, tol: float = 0.04) -> pd.DataFrame:
    """Same band in different lanes -> same ``band_group`` (1 = top of the gel).

    Bands are compared by size when the gel is calibrated (within ``tol``, relative) and by
    relative migration (Rf, within ``tol`` / 2) otherwise. A group holds at most one band per
    lane. Marker lanes are not grouped."""
    df = df.copy()
    df["band_group"] = pd.array([pd.NA] * len(df), dtype="Int64")
    df["group_size"] = np.nan
    if not len(df):
        return df
    use = df[~df["is_marker"]].copy()
    by_size = size_col is not None and size_col in use and use[size_col].notna().all() and len(use)
    if by_size:
        key = -np.log10(use[size_col].astype(float))
        tol_key = np.log10(1 + tol)
    else:
        key = use["rf"].astype(float)
        tol_key = tol / 2
    order = np.argsort(key.to_numpy(), kind="stable")
    idx = use.index.to_numpy()[order]
    keys = key.to_numpy()[order]
    lanes = use["lane"].to_numpy()[order]
    groups: list[list[int]] = []
    cur, cur_keys, cur_lanes = [], [], set()
    for i, k, ln in zip(idx, keys, lanes):
        if cur and (abs(k - np.mean(cur_keys)) > tol_key or ln in cur_lanes):
            groups.append(cur)
            cur, cur_keys, cur_lanes = [], [], set()
        cur.append(i)
        cur_keys.append(k)
        cur_lanes.add(ln)
    if cur:
        groups.append(cur)
    for g, members in enumerate(groups, start=1):
        df.loc[members, "band_group"] = g
        if by_size:
            df.loc[members, "group_size"] = round(float(np.median(df.loc[members, size_col])), 1)
    return df


def normalise(df: pd.DataFrame, ga, size_col: str | None) -> tuple[pd.DataFrame, list[str]]:
    """Western-blot style normalisation (Params.norm_method):

    * ``loading_control``: divide each band by the lane's loading-control band (the strongest
      band within ``lc_tol`` of ``lc_size``, e.g. 42 kDa for actin), optionally read from another
      channel (``lc_channel``, two-colour blots).
    * ``total_lane``: divide by the lane's total signal (total-protein / stain-free / Ponceau),
      optionally from another channel or an added image of the total-protein stain.

    ``fold_change`` = normalised value / the same band group's value in ``reference_lane``."""
    p = ga.params
    warn: list[str] = []
    df = df.copy()
    if p.norm_method in (None, "", "none") or not len(df):
        return df, warn
    names = ga.channel_names()
    lc_col = tgt_col = "volume"
    for attr, label in (("lc_channel", "loading-control"), ("target_channel", "target")):
        ch = getattr(p, attr)
        if ch is not None:
            if not 0 <= ch < len(names) or len(names) < 2:
                return df, [f"{label} channel {ch + 1} does not exist"]
            if attr == "lc_channel":
                lc_col = channel_column(names[ch])
            else:
                tgt_col = channel_column(names[ch])
    lc = {}
    if p.norm_method == "loading_control":
        if not p.lc_size or size_col is None or df[size_col].isna().all():
            return df, ["Normalisation: set the loading-control size and calibrate with a marker first."]
        if lc_col not in df:
            return df, [f"no volumes for channel column {lc_col}"]
        df["is_loading_control"] = False
        for lane, sub in df[~df["is_marker"]].groupby("lane"):
            cand = sub[(sub[size_col] - p.lc_size).abs() <= p.lc_tol * p.lc_size]
            if len(cand):
                k = cand[lc_col].idxmax()
                lc[lane] = float(df.loc[k, lc_col])
                df.loc[k, "is_loading_control"] = True
            else:
                warn.append(f"lane {lane}: no loading-control band near {p.lc_size:g}")
    elif p.norm_method == "total_lane":
        for i in range(len(ga.lanes)):
            if (i + 1) in p.marker_lanes:
                continue
            pr = ga.lane_profile(i, channel=p.lc_channel if p.lc_channel is not None else "auto")
            lc[i + 1] = float(np.clip(pr["corr"], 0, None).sum())
    else:
        return df, [f"unknown normalisation {p.norm_method!r}"]
    df["lc_signal"] = df["lane"].map(lc)
    df["normalised"] = df[tgt_col] / df["lc_signal"]
    df.loc[df["is_marker"], "normalised"] = np.nan
    if p.reference_lane:
        ref = df[(df["lane"] == p.reference_lane) & df["band_group"].notna()]
        ref_val = dict(zip(ref["band_group"].astype(int), ref["normalised"]))
        if not ref_val:
            warn.append(f"reference lane {p.reference_lane} has no bands to compare with")
        df["fold_change"] = [v / ref_val[int(g)] if pd.notna(g) and int(g) in ref_val and ref_val[int(g)] else np.nan
                             for v, g in zip(df["normalised"], df["band_group"])]
    return df, warn


def group_matrix(df: pd.DataFrame, value: str = "volume") -> pd.DataFrame:
    """Band groups x lanes table of ``value`` (e.g. volume, normalised, fold_change)."""
    if "band_group" not in df or value not in df or df["band_group"].isna().all():
        return pd.DataFrame()
    sub = df[df["band_group"].notna()]
    m = sub.pivot_table(index="band_group", columns="lane_name", values=value, aggfunc="first")
    order = list(dict.fromkeys(sub.sort_values("lane")["lane_name"]))
    m = m.reindex(columns=order)
    sizes = sub.groupby("band_group")["group_size"].first()
    m.insert(0, "size", sizes)
    return m.reset_index()


# RNA: typical eukaryotic rRNA sizes (nt); 28S 4.7-5.1 kb in mammals, 3.4 kb in plants (25S)
RNA_28S = (3200.0, 6000.0)
RNA_18S = (1500.0, 2300.0)


def rna_integrity(ga, df: pd.DataFrame, size_col: str | None) -> pd.DataFrame:
    """Per lane: 28S/18S volume ratio and the share of signal below the 18S band (degradation).

    With a calibrated gel the rRNA bands are found by size; otherwise the two strongest bands are
    taken (upper = 28S, lower = 18S). Intact mammalian total RNA has a ratio around 1.8-2.2 and
    little signal below 18S; degraded RNA shows a lower ratio and a smear towards small sizes."""
    rows = []
    for i, lane in enumerate(ga.lanes):
        if (i + 1) in ga.params.marker_lanes:
            continue
        sub = df[df["lane"] == i + 1]
        r = {"lane": i + 1, "lane_name": lane.name, "rna_28s_18s": np.nan, "rna_below_18s_pct": np.nan,
             "rna_note": ""}
        b28 = b18 = None
        if size_col and len(sub) and sub[size_col].notna().any():
            c28 = sub[sub[size_col].between(*RNA_28S)]
            c18 = sub[sub[size_col].between(*RNA_18S)]
            if len(c28) and len(c18):
                b28, b18 = c28.loc[c28["volume"].idxmax()], c18.loc[c18["volume"].idxmax()]
        if b28 is None and len(sub) >= 2:
            two = sub.nlargest(2, "volume").sort_values("y_px")
            b28, b18 = two.iloc[0], two.iloc[1]
            r["rna_note"] = "rRNA bands taken as the two strongest bands"
        if b28 is None:
            r["rna_note"] = "fewer than two bands"
            rows.append(r)
            continue
        r["rna_28s_18s"] = round(float(b28["volume"] / b18["volume"]), 2) if b18["volume"] > 0 else np.nan
        y, sig = ga.smear_profile(i)          # background beside the lane: keeps the smear
        top = y >= b28["y0_px"]
        total = float(sig[top].sum())
        below = float(sig[y > b18["y1_px"]].sum())
        r["rna_below_18s_pct"] = round(100 * below / total, 1) if total > 0 else np.nan
        rows.append(r)
    return pd.DataFrame(rows)

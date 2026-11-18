"""Comparing lanes: channels / added images, band groups, Western-blot normalisation, RNA integrity.

Simulated gels with known truth:
* Western blot: target band (55 kDa) with amounts 0.5-4x, loading control (42 kDa) with uneven
  loading; the chemiluminescence image shows no marker - it is in a second (white-light) image.
* two-colour blot: target in page 1, loading control in page 2 of one TIFF.
* RNA: intact lanes (28S:18S = 2) and degraded lanes (ratio 1 plus a smear below 18S).
"""
import numpy as np
import pytest
import tifffile

from autogel.analysis import GelAnalysis, Params

from .synth import make_gel

PROT = [250, 150, 100, 75, 50, 37, 25, 20, 15, 10]       # kDa, top -> bottom


def y_kda(s, y_top=200, y_bottom=1100):
    return y_top + (np.log10(PROT[0]) - np.log10(s)) / (np.log10(PROT[0]) - np.log10(PROT[-1])) * (y_bottom - y_top)


AMOUNT = {1: 1.0, 2: 1.0, 3: 2.0, 4: 2.0, 5: 4.0, 6: 4.0, 7: 0.5, 8: 0.5, 9: 1.0}   # lane index -> target
LOAD = {1: 1.0, 2: 0.8, 3: 1.2, 4: 1.0, 5: 0.9, 6: 1.1, 7: 1.0, 8: 1.3, 9: 0.7}     # protein loaded


def _western(tmp_path):
    lanes = {k: [(y_kda(55), 3000 * AMOUNT[k] * LOAD[k], 5), (y_kda(42), 6000 * LOAD[k], 5)] for k in AMOUNT}
    chemi = make_gel(str(tmp_path / "chemi.tif"), lanes=lanes, bg=800, gradient=200, noise=30)
    marker = make_gel(str(tmp_path / "marker.tif"), lanes={0: [(y_kda(s), 9000, 4) for s in PROT]},
                      invert=True, seed=3, noise=30)        # dark prestained bands on a light membrane
    return chemi, marker


def test_marker_from_separate_image_and_loading_control(tmp_path):
    chemi, marker = _western(tmp_path)
    ga = GelAnalysis(chemi, Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    assert sum(b.lane == 0 for b in ga.bands) == 0            # no marker visible in the chemi image
    k = ga.add_image(marker, name="white light")
    assert ga.channel_names() == ["image", "white light"] and k == 1
    assert ga.params.extra_images[0]["invert"] is True       # light background recognised
    ga.set_marker([1], sizes=PROT, unit="kDa")
    ga.set_marker_channel(1)
    assert sum(b.lane == 0 for b in ga.bands) == len(PROT)   # read from the marker image
    ga.params.norm_method = "loading_control"
    ga.params.lc_size = 42
    ga.params.reference_lane = 2
    df = ga.quantify()
    tgt = df[(~df.is_marker) & (df.size_kDa.between(50, 60))]
    assert len(tgt) == 9 and (tgt.size_kDa.sub(55).abs() < 55 * 0.03).all()
    assert tgt.band_group.nunique() == 1                      # same band in all lanes
    lc = df[df.is_loading_control]
    assert len(lc) == 9 and (lc.size_kDa.sub(42).abs() < 42 * 0.03).all()
    fc = dict(zip(tgt.lane, tgt.fold_change))
    for k, amount in AMOUNT.items():
        assert fc[k + 1] == pytest.approx(amount / AMOUNT[1], rel=0.06), k
    # session replay re-attaches the marker image and keeps the marker channel
    s = tmp_path / "s.json"
    ga.save_session(str(s))
    gb = GelAnalysis.load_session(str(s))
    assert gb.channel_names() == ["image", "white light"] and gb.params.marker_channel == 1
    np.testing.assert_allclose(gb.quantify()["fold_change"], df["fold_change"], rtol=1e-6)


def test_added_image_must_match(tmp_path):
    chemi, _ = _western(tmp_path)
    other = make_gel(str(tmp_path / "small.tif"), H=600, W=500, lanes={})
    ga = GelAnalysis(chemi, Params(n_lanes=10))
    with pytest.raises(ValueError, match="same camera position"):
        ga.add_image(other)
    with pytest.raises(ValueError):
        ga.set_channel(3)
    with pytest.raises(ValueError):
        ga.set_marker_channel(2)


def test_total_lane_normalisation(tmp_path):
    chemi, _ = _western(tmp_path)
    ga = GelAnalysis(chemi, Params(n_lanes=10, norm_method="total_lane", reference_lane=2))
    ga.detect_lanes()
    ga.detect_bands()
    df = ga.quantify()
    assert "normalised" in df and "fold_change" in df
    assert df[df.lane == 2]["fold_change"].dropna().round(6).eq(1).all()


def test_two_colour_blot(tmp_path):
    """Target (page 1) and loading control (page 2) in one two-page TIFF."""
    t = make_gel(str(tmp_path / "t.tif"), lanes={k: [(y_kda(55), 3000 * AMOUNT[k] * LOAD[k], 5)] for k in AMOUNT},
                 noise=30)
    c = make_gel(str(tmp_path / "c.tif"), lanes={k: [(y_kda(42), 6000 * LOAD[k], 5)] for k in AMOUNT}, noise=30,
                 seed=5)
    m = make_gel(str(tmp_path / "m.tif"), lanes={0: [(y_kda(s), 9000, 4) for s in PROT]}, noise=30, seed=7)
    two = str(tmp_path / "two.tif")
    tifffile.imwrite(two, np.stack([tifffile.imread(t), tifffile.imread(c)]))
    ga = GelAnalysis(two, Params(n_lanes=10))
    assert ga.channel_names() == ["page 1", "page 2"]
    assert any("Multi-page" in w for w in ga.gel.warnings)
    ga.add_image(m, name="marker")
    ga.set_channel(-1)                                         # detect bands of both colours
    ga.detect_lanes()
    ga.detect_bands()
    ga.set_marker([1], sizes=PROT, unit="kDa")
    ga.set_marker_channel(2)
    p = ga.params
    p.norm_method, p.lc_size, p.lc_channel, p.target_channel, p.reference_lane = "loading_control", 42, 1, 0, 2
    df = ga.quantify()
    assert {"volume_page_1", "volume_page_2", "volume_marker"} <= set(df.columns)
    tgt = df[(~df.is_marker) & df.size_kDa.between(50, 60)]
    fc = dict(zip(tgt.lane, tgt.fold_change))
    for k, amount in AMOUNT.items():
        assert fc[k + 1] == pytest.approx(amount / AMOUNT[1], rel=0.06), k


RNA_LAD = [9000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500]


def y_nt(s, y_top=200, y_bottom=1100):
    return y_top + (np.log10(RNA_LAD[0]) - np.log10(s)) / (np.log10(RNA_LAD[0]) - np.log10(RNA_LAD[-1])) \
        * (y_bottom - y_top)


def test_rna_integrity(tmp_path):
    lanes = {0: [(y_nt(s), 6000, 4) for s in RNA_LAD]}
    intact = (1, 2, 3, 4)
    degraded = (5, 6, 7, 8)
    for k in intact:
        lanes[k] = [(y_nt(5000), 8000, 5), (y_nt(1900), 4000 * 1.05, 5)]    # areas 2:1 (18S slightly broader)
    for k in degraded:
        smear = [(y, 700, 12) for y in np.arange(y_nt(1700), y_nt(500), 15)]
        lanes[k] = [(y_nt(5000), 4000, 5), (y_nt(1900), 4000, 5)] + smear
    p = make_gel(str(tmp_path / "rna.tif"), lanes=lanes, noise=30)
    ga = GelAnalysis(p, Params(n_lanes=10, rna_integrity=True, min_coverage=0.2))
    ga.detect_lanes()
    ga.detect_bands()
    ga.set_marker([1], sizes=RNA_LAD, unit="nt")
    ga.quantify()
    lt = ga.lanes_table().set_index("lane")
    for k in intact:
        assert lt.loc[k + 1, "rna_28s_18s"] == pytest.approx(2.0, rel=0.15)
    for k in degraded:
        assert lt.loc[k + 1, "rna_28s_18s"] < 1.4
        assert lt.loc[k + 1, "rna_below_18s_pct"] > 2 * lt.loc[intact[0] + 1, "rna_below_18s_pct"] + 5


def test_band_groups_without_calibration(tmp_path):
    """Without a marker, bands are matched by relative migration."""
    lanes = {k: [(400, 5000, 4), (700, 3000, 4)] + ([(900, 2000, 4)] if k % 2 else []) for k in range(1, 9)}
    ga = GelAnalysis(make_gel(str(tmp_path / "g.tif"), lanes=lanes), Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    df = ga.quantify()
    sub = df[df.lane.between(2, 9)]
    assert sub.groupby("band_group").lane.nunique().tolist() == [8, 8, 4]
    from autogel.compare import group_matrix
    m = group_matrix(df)
    assert list(m.columns[:2]) == ["band_group", "size"] and len(m) == 3

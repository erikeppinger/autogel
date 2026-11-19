import os

import numpy as np
import pytest
import tifffile

from autogel.analysis import Band, GelAnalysis, Lane, Params
from autogel.export import export_all

from .synth import ladder, make_gel, y_for_size

SIZES = [10000, 6000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 100]


def standard_lanes():
    lad, _ = ladder(sizes=SIZES)
    lanes = {0: lad, 9: lad}
    for k in range(1, 9):
        lanes[k] = [(y_for_size(2500), 6000, 5), (y_for_size(800), 3000, 5)]
    return lanes


@pytest.fixture
def gel(tmp_path):
    return make_gel(str(tmp_path / "gel.tif"), lanes=standard_lanes())


def run(path, **kw):
    ga = GelAnalysis(path, Params(**kw))
    ga.detect_lanes()
    ga.detect_bands()
    return ga


# ------------------------------------------------------------------ lanes
def test_auto_lane_count(gel):
    ga = run(gel)
    assert len(ga.lanes) == 10
    xs = [l.x for l in ga.lanes]
    assert np.allclose(xs, [140 + 80 * k for k in range(10)], atol=4)


def test_given_lane_count(gel):
    assert len(run(gel, n_lanes=10).lanes) == 10


def test_empty_lane_is_filled(tmp_path):
    lanes = standard_lanes()
    lanes[4] = []
    ga = run(make_gel(str(tmp_path / "g.tif"), lanes=lanes), n_lanes=10)
    assert len(ga.lanes) == 10
    assert abs(ga.lanes[4].x - (140 + 4 * 80)) < 6
    assert not any(b.lane == 4 for b in ga.bands)


def test_rotated_input(tmp_path):
    ga = run(make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), rot90=1), n_lanes=10)
    assert len(ga.lanes) == 10
    assert sum(b.lane == 0 for b in ga.bands) == len(SIZES)


def test_upside_down_input(tmp_path):
    ga = run(make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), rot90=2), n_lanes=10)
    # ladder still reads large -> small from the top
    lane0 = sorted(b.y for b in ga.bands if b.lane == 0)
    assert len(lane0) == len(SIZES)
    ga.set_marker([1], sizes=SIZES)
    df = ga.quantify()
    assert df[df.lane == 2]["size_bp"].iloc[0] == pytest.approx(2500, rel=0.05)


def test_inverted_stain(tmp_path):
    ga = run(make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), invert=True), n_lanes=10)
    assert ga.gel.orientation["invert"] is True
    assert len(ga.lanes) == 10
    assert sum(b.lane == 0 for b in ga.bands) == len(SIZES)


# ------------------------------------------------------------------ bands
def test_ladder_band_count_and_no_false_bands(gel):
    ga = run(gel, n_lanes=10)
    counts = [sum(b.lane == i for b in ga.bands) for i in range(10)]
    assert counts[0] == counts[9] == len(SIZES)
    assert counts[1:9] == [2] * 8


def test_relative_volume(gel):
    ga = run(gel, n_lanes=10)
    df = ga.quantify()
    lane = df[df.lane == 3].sort_values("y_px")
    ratio = lane.volume.iloc[0] / lane.volume.iloc[1]
    assert ratio == pytest.approx(2.0, rel=0.05)    # 6000 vs 3000 amplitude, same sigma


def test_dust_is_rejected(tmp_path):
    p = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), dust=[(900, 140 + 3 * 80, 20000)])
    ga = run(p, n_lanes=10)
    assert sum(b.lane == 3 for b in ga.bands) == 2


def test_saturation_flag(tmp_path):
    lanes = standard_lanes()
    lanes[5] = [(y_for_size(2500), 70000, 5)]
    p = make_gel(str(tmp_path / "g.tif"), lanes=lanes, clip=60000)
    ga = run(p, n_lanes=10)
    df = ga.quantify()
    assert df[df.lane == 6].saturated.all()
    assert not df[df.lane == 4].saturated.any()
    assert ga.gel.saturation_level == 60000


# ------------------------------------------------------------ calibration
def test_size_calibration(gel):
    ga = run(gel, n_lanes=10)
    ga.set_marker([1, 10], sizes=SIZES)
    df = ga.quantify()
    s = df[(df.lane == 5)].sort_values("y_px")["size_bp"].tolist()
    assert s[0] == pytest.approx(2500, rel=0.03)
    assert s[1] == pytest.approx(800, rel=0.03)
    assert not df[df.lane == 5].size_extrapolated.any()


def test_marker_count_mismatch_warns(gel):
    ga = run(gel, n_lanes=10)
    ga.set_marker([1], sizes=SIZES + [50])          # 50 bp is not on the gel
    df = ga.quantify()
    warns = [w for c in ga.calibrations for w in c.warnings]
    assert any("assumed missing: 50" in w for w in warns)
    # the best-fit assignment keeps the sizes right despite the extra ladder entry
    assert df[df.lane == 2]["size_bp"].iloc[0] == pytest.approx(2500, rel=0.03)


def _neb_lane():
    """NEB 1 kb Plus as on a real gel: rows from a curved migration law, volumes ~ masses with a
    smooth size dependence (as measured on the GelGenie ladder gels)."""
    from autogel.markers import PRESETS
    e = PRESETS["NEB Quick-Load Purple 1 kb Plus (N0550)"]
    ls = np.log10(e["sizes"])
    x = (ls.max() - ls) / (ls.max() - ls.min())
    ys = 240 + 560 * (0.7 * x + 0.3 * x ** 2)
    vols = np.asarray(e["masses"], float) * 100 * (1 - 0.4 * x)
    return e, list(ys), list(vols)


def test_marker_wells_and_dye_front_beyond_three_extras():
    # 2 bands at the wells + 2 at the dye front: the old search gave up above 3 extras and
    # paired top -> bottom, so every size was wrong (24_NEB lane 6 of the GelGenie gels)
    from autogel.markers import assign_marker
    e, ys, vols = _neb_lane()
    yy, vv = [40.0, 55.0] + ys + [930.0, 955.0], [3000.0, 2000.0] + vols + [4000.0, 3000.0]
    y_a, s_a, _, warn = assign_marker(yy, e["sizes"], e["masses"], vv)
    assert dict(zip(y_a, s_a)) == dict(zip(ys, e["sizes"]))
    assert any("Ignored the band" in w for w in warn)


def test_marker_merged_pair_plus_well_band():
    # 1000 bp swallows 900 bp (merged, at the bright band's place) and a well band is detected:
    # same count as the ladder, so the old code paired directly and shifted the top
    from autogel.markers import assign_marker
    e, ys, vols = _neb_lane()
    k = e["sizes"].index(1000)
    yy = [90.0] + ys[:k + 1] + ys[k + 2:]
    vv = [5000.0] + vols[:k] + [vols[k] + vols[k + 1]] + vols[k + 2:]
    y_a, s_a, _, _ = assign_marker(yy, e["sizes"], e["masses"], vv)
    got = dict(zip(y_a, s_a))
    assert 90.0 not in got                                  # the well band is not a ladder band
    truth = dict(zip(ys, e["sizes"]))
    wrong = [y for y, s in got.items() if y != ys[k] and truth[y] != s]
    assert not wrong
    assert got.get(ys[k]) in (1000, 900, None)


def _shoulder_gel(tmp_path, sig=9):
    """NEB ladder whose 1000 and 900 bp bands are wide, so 900 is only a shoulder of 1000."""
    from autogel.markers import PRESETS
    e = PRESETS["NEB Quick-Load Purple 1 kb Plus (N0550)"]
    bands, ymap = ladder(sizes=tuple(e["sizes"]), amp=1, sig=4)
    w = [sig if s in (1000, 900) else 4 for s in e["sizes"]]
    lanes = [(y, 480 * m / s_, s_) for (y, _, _), m, s_ in zip(bands, e["masses"], w)]
    return make_gel(str(tmp_path / "shoulder.tif"), n_lanes=3, lanes={0: lanes, 1: lanes, 2: lanes}), ymap


def test_ladder_guided_split_finds_shoulder_band(tmp_path):
    path, ymap = _shoulder_gel(tmp_path)
    neb = "NEB Quick-Load Purple 1 kb Plus (N0550)"
    found = {}
    for guided in (False, True):
        ga = GelAnalysis(path, Params(n_lanes=3, ladder_guided=guided))
        ga.detect_lanes()
        ga.detect_bands()
        ga.set_marker([1, 2], preset=neb)
        ga.quantify()
        c = ga.calibrations[0]
        found[guided] = dict(zip(c.sizes, c.ys))
    assert 900 not in found[False]                          # without guidance: a shoulder only
    assert found[True][900] == pytest.approx(ymap[900], abs=5)
    # the split part is marked and not used as an ng standard
    split = [b for b in ga.bands if b.source == "ladder"]
    assert split and all(abs(b.y - ymap[900]) < 6 for b in split)
    assert any(e["step"] == "refine_marker_bands" for e in ga.log)


def test_ladder_split_respects_user_deletion_and_sessions(tmp_path):
    path, _ = _shoulder_gel(tmp_path)
    ga = GelAnalysis(path, Params(n_lanes=3))
    ga.detect_lanes()
    ga.detect_bands()
    ga.set_marker([1], preset="NEB Quick-Load Purple 1 kb Plus (N0550)")
    ga.quantify()
    assert any(b.source == "ladder" for b in ga.bands)
    ga.set_bands([b for b in ga.bands if b.source != "ladder"], reason="user deleted the split band")
    ga.quantify()                                           # must not come back
    assert not any(b.source == "ladder" for b in ga.bands)
    sp = str(tmp_path / "s.json")
    ga.save_session(sp)
    ga2 = GelAnalysis.load_session(sp)
    ga2.quantify()
    assert not any(b.source == "ladder" for b in ga2.bands)
    ga2.detect_bands(lanes=[0])                             # a new detection refines again
    ga2.quantify()
    assert any(b.source == "ladder" for b in ga2.bands)


def test_marker_with_dust_band_is_matched(gel):
    ga = run(gel, n_lanes=10)
    ga.set_marker([1], sizes=SIZES)
    lane0 = [b for b in ga.bands if b.lane == 0]
    fake = ga.band_for(0, lane0[0].y - 60, lane0[0].y - 50)   # spurious band above 10 kb
    ga.set_bands(ga.bands + [fake])
    df = ga.quantify()
    assert any("Ignored the band" in w for c in ga.calibrations for w in c.warnings)
    assert df[df.lane == 2]["size_bp"].iloc[0] == pytest.approx(2500, rel=0.03)


def test_mass_curve_and_load_scaling(gel):
    ga = run(gel, n_lanes=10)
    masses = [10.0] * len(SIZES)       # all ladder bands equal -> equal volumes
    ga.set_marker([1, 10], sizes=SIZES, masses=masses)
    df1 = ga.quantify()
    m1 = df1[df1.lane == 3]["est_mass_ng"].iloc[0]
    # volume ~ amplitude x sigma: sample 6000 x 5 vs ladder 8000 x 4 per 10 ng
    assert m1 == pytest.approx(10 * 6000 * 5 / (8000 * 4), rel=0.08)
    ga.set_marker([1, 10], sizes=SIZES, masses=masses, load=2, ref_load=1)   # twice the ladder loaded
    df2 = ga.quantify()
    assert df2[df2.lane == 3]["est_mass_ng"].iloc[0] == pytest.approx(2 * m1, rel=2e-3)
    ga.params.sample_loads_ul = {"3": 5}
    df3 = ga.quantify()
    assert df3[df3.lane == 3]["conc_ng_per_ul"].iloc[0] == pytest.approx(2 * m1 / 5, rel=1e-3)


# ------------------------------------------------------------- user edits
def test_manual_band_and_lane_edits(gel):
    ga = run(gel, n_lanes=10)
    n = len(ga.bands)
    b = ga.band_at(4, y_for_size(2500) + 1)
    assert abs(b.y - y_for_size(2500)) < 2
    ga.set_bands(ga.bands + [ga.band_for(4, 1000, 1010)])
    assert len(ga.bands) == n + 1
    assert ga.log[-1]["step"] == "set_bands"
    lanes = [Lane(x=l.x + 3, width=l.width, y0=l.y0, y1=l.y1) for l in ga.lanes]
    ga.set_lanes(lanes)
    assert len(ga.bands) == n + 1          # bands follow their (moved) lane
    ga.set_lanes(lanes[:-1])               # delete last lane
    assert not any(b.lane == 9 for b in ga.bands)


def test_session_roundtrip_reproduces_numbers(gel, tmp_path):
    ga = run(gel, n_lanes=10)
    ga.set_bands(ga.bands + [Band(lane=4, y=1005, y0=1000, y1=1010, source="user")])
    ga.set_marker([1, 10], sizes=SIZES)
    df = ga.quantify()
    s = str(tmp_path / "s.json")
    ga.save_session(s)
    ga2 = GelAnalysis.load_session(s)
    df2 = ga2.quantify()
    assert df.volume.tolist() == df2.volume.tolist()
    assert df.size_bp.tolist() == pytest.approx(df2.size_bp.tolist(), nan_ok=True)
    assert any(e["step"] == "load_session" for e in ga2.log)


def test_export_files(gel, tmp_path):
    ga = run(gel, n_lanes=10)
    ga.set_marker([1, 10], sizes=SIZES, masses=[10.0] * len(SIZES))
    out = export_all(ga, str(tmp_path / "out"))
    names = set(os.listdir(out))
    for suffix in ("_bands.csv", "_lanes.csv", "_profiles.csv", "_profiles.pdf", "_annotated.png",
                   "_log.txt", "_session.json", "_autogel.xlsx"):
        assert "gel" + suffix in names


# ---------------------------------------------------------- break tests
def test_blank_image_does_not_crash(tmp_path):
    p = str(tmp_path / "blank.tif")
    tifffile.imwrite(p, np.full((400, 500), 1000, np.uint16))
    ga = run(p)
    assert ga.lanes == [] or isinstance(ga.lanes, list)
    df = ga.quantify()
    assert len(df) == 0


def test_noise_only_image(tmp_path):
    p = str(tmp_path / "noise.tif")
    tifffile.imwrite(p, np.random.default_rng(1).normal(1000, 50, (400, 500)).clip(0).astype(np.uint16))
    ga = run(p)
    ga.quantify()


def test_single_lane(tmp_path):
    lad, _ = ladder(sizes=SIZES)
    p = make_gel(str(tmp_path / "g.tif"), n_lanes=1, W=300, lanes={0: lad})
    ga = run(p, n_lanes=1)
    assert len(ga.lanes) == 1
    assert sum(b.lane == 0 for b in ga.bands) == len(SIZES)


def test_more_lanes_requested_than_visible(gel):
    ga = run(gel, n_lanes=12)
    assert len(ga.lanes) == 12


def test_8bit_and_rgb_inputs(tmp_path):
    from skimage import io
    p16 = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes())
    a = tifffile.imread(p16)
    a8 = (a / 257).astype(np.uint8)
    p8 = str(tmp_path / "g8.png")
    io.imsave(p8, np.dstack([a8, a8, a8]), check_contrast=False)
    ga = run(p8, n_lanes=10)
    assert len(ga.lanes) == 10
    assert any("8-bit" in w for w in ga.gel.warnings)


def test_band_outside_profile_range(gel):
    ga = run(gel, n_lanes=10)
    b = ga.band_for(2, -50, -10)      # nonsense limits from a careless edit
    ga.set_bands(ga.bands + [b])
    ga.quantify()


def test_bad_marker_lane_numbers(gel):
    ga = run(gel, n_lanes=10)
    ga.set_marker([0, 99], sizes=SIZES)
    df = ga.quantify()
    assert df["size_bp"].isna().all()


def test_degenerate_user_lanes(gel):
    """Lanes dragged outside the image, zero height or zero width must not crash."""
    ga = run(gel, n_lanes=10)
    H, W = ga.data.shape
    lanes = list(ga.lanes) + [
        Lane(x=W + 50, width=40, y0=100, y1=900),        # completely outside
        Lane(x=W - 5, width=40, y0=100, y1=900),         # half outside
        Lane(x=500, width=40, y0=600, y1=600),           # zero height
        Lane(x=520, width=0.1, y0=100, y1=900),          # zero width
        Lane(x=300, width=40, y0=900, y1=100),           # upside-down limits
    ]
    ga.set_lanes(lanes)
    ga.detect_bands()
    ga.quantify()
    ga.lanes_table()
    for i in range(len(ga.lanes)):
        ga.band_at(i, 500)


def test_session_with_missing_image(gel, tmp_path):
    ga = run(gel, n_lanes=10)
    s = str(tmp_path / "s.json")
    ga.save_session(s)
    os.remove(gel)
    with pytest.raises(FileNotFoundError):
        GelAnalysis.load_session(s)


def test_user_ladder_library(tmp_path, monkeypatch):
    monkeypatch.setenv("AUTOGEL_HOME", str(tmp_path))
    from autogel import markers
    markers.save_user_ladder("my ladder", [100, 500, 1000], "bp", [10, 20, 30])
    assert markers.PRESETS["my ladder"]["sizes"] == [1000, 500, 100]
    assert markers.PRESETS["my ladder"]["user_defined"] is True
    with pytest.raises(ValueError):
        markers.save_user_ladder("bad", [100, 200], "bp", [1.0])
    monkeypatch.delenv("AUTOGEL_HOME")
    markers.PRESETS.clear()
    markers.PRESETS.update(markers.load_library())


def test_label_format_and_declutter(gel):
    from autogel.labels import declutter, fmt_size, plan_labels
    assert [fmt_size(s) for s in (10002, 8001, 1517, 1000, 900)] == ["10 kb", "8.0 kb", "1.5 kb", "1.0 kb", "900"]
    assert fmt_size(37, "kDa") == "37"
    ga = run(gel, n_lanes=10)
    ga.set_marker([1, 10], sizes=SIZES)
    df = ga.quantify()
    labs = plan_labels(ga, df, {0, 9, 4})
    sides = {(l.lane, l.side) for l in labs}
    assert (0, "left") in sides and (9, "right") in sides and (4, "right") in sides
    for l in labs:                           # never on the band: outside the lane edges
        lane = ga.lanes[l.lane]
        assert l.x < lane.x0 or l.x > lane.x1
    keep_far = declutter(labs, 1.0)
    keep_near = declutter(labs, 400.0)
    assert all(keep_far) and sum(keep_near) < len(labs)
    # each lane side keeps at least its highest-priority label
    assert {(l.lane, l.side) for l, k in zip(labs, keep_near) if k} == sides


def test_size_range_limits_calibration_and_flags(gel):
    ga = run(gel, n_lanes=10)
    ga.set_marker([1, 10], sizes=SIZES)
    df_all = ga.quantify()
    ga.set_size_range(500, 3000)
    df = ga.quantify()
    assert all(500 * 0.98 <= s <= 3000 * 1.02 for s in ga.calibrations[0].sizes)
    s2500 = df[(df.lane == 3)].sort_values("y_px").iloc[0]
    assert s2500["size_in_range"] and s2500["size_bp"] == pytest.approx(2500, rel=0.03)
    assert df.size_in_range.sum() < df_all.size_in_range.sum()
    ga.set_size_range(9000, 9500)          # < 2 ladder bands inside -> warning, range ignored
    ga.quantify()
    assert any("fewer than 2" in w for c in ga.calibrations for w in c.warnings)


@pytest.mark.parametrize("rot", [0, 1])
def test_orientation_overrules_smear_fooled_band_shape(tmp_path, monkeypatch, rot):
    """Heavily smeared lanes make band blobs look tall (band-shape vote says 'rotate'); the
    well row and the lane comb must keep the decision right, as on the real benchmark gels."""
    import autogel.imageio as io
    path = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), rot90=rot)
    a = np.asarray(tifffile.imread(path))
    truth = io.detect_orientation(a)["rot90_k"] % 2
    assert truth == rot % 2
    monkeypatch.setattr(io, "lanes_vertical_vote", lambda s: 0.3)   # smear-fooled, says "turn 90°"
    o = io.detect_orientation(a)
    assert o["rot90_k"] % 2 == rot % 2
    assert o["vertical_cues"]["band_shape"] == pytest.approx(-0.2)


def test_uncertain_orientation_warns(tmp_path, monkeypatch):
    import autogel.imageio as io
    path = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes())
    real = io.detect_orientation
    monkeypatch.setattr(io, "detect_orientation", lambda a, inv=None: {**real(a, inv), "vertical_score": 0.02})
    assert any("Lane direction uncertain" in w for w in io.load_gel(path).warnings)
    # an orientation given by the user (Rotate 90°, saved session) is not questioned
    o = real(np.asarray(tifffile.imread(path)))
    assert not any("uncertain" in w for w in io.load_gel(path, orientation=o).warnings)


def test_skew_reestimated_after_manual_rotation(tmp_path):
    """Rotate 90° in the GUI sets angle=None: the skew is measured again in the new frame,
    not carried over from the wrong one."""
    from autogel.imageio import load_gel
    path = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), rot90=1)
    auto = load_gel(path).orientation
    g = load_gel(path, orientation={**auto, "rot90_k": (auto["rot90_k"] + 2) % 4, "angle": None})
    assert isinstance(g.orientation["angle"], float) and abs(g.orientation["angle"]) <= 3.0
    assert not any("uncertain" in w for w in g.warnings)


def test_size_dependent_mass_curve(tmp_path):
    """Signal per ng falls for small fragments (as on the real gels: ~20 % at 100 bp). The size
    model learns that from the ladder; a single slope over-/under-estimates by several fold."""
    def resp(size):                       # relative signal per ng, peak at ~2 kb
        return np.exp(-0.5 * ((np.log10(size) - 3.3) / 0.55) ** 2)
    masses = [10.0, 20, 15, 30, 10, 25, 40, 10, 20, 15, 30]
    lad = [(y_for_size(s), 3000 * m / 10 * resp(s), 4) for s, m in zip(SIZES, masses)]
    true = {2: (250, 20.0), 3: (2500, 20.0), 4: (6000, 8.0)}            # lane: (size, ng)
    lanes = {0: lad, 9: lad}
    for k, (size, m) in true.items():
        lanes[k] = [(y_for_size(size), 3000 * m / 10 * resp(size), 4)]
    ga = run(make_gel(str(tmp_path / "g.tif"), lanes=lanes, noise=20), n_lanes=10)
    ga.set_marker([1, 10], sizes=SIZES, masses=masses)
    df = ga.quantify()
    assert ga.mass_fit["model"] == "size" and ga.mass_fit["typical_error"] < 0.15
    for k, (size, m) in true.items():
        row = df[df.lane == k + 1].iloc[0]
        assert row["mass_size_corrected"]
        assert row["est_mass_ng"] == pytest.approx(m, rel=0.15), (size, row["est_mass_ng"])
    ga.params.mass_model = "linear"
    lin = ga.quantify()
    assert ga.mass_fit["model"] == "linear" and not lin["mass_size_corrected"].any()
    small = lin[lin.lane == 3].iloc[0]["est_mass_ng"]
    assert small < 0.5 * 20.0                                           # one slope: far too low at 250 bp


def _doublet_gel(tmp_path, sep, ratio, name="d.tif"):
    lanes = {k: [(500, 6000, 4.0), (500 + sep, 6000 / ratio, 4.0)] for k in range(1, 9)}
    lanes[5] = [(500, 6000, 6.0)]                         # a broad single band must stay one band
    return make_gel(str(tmp_path / name), lanes=lanes, noise=30)


def _ratios(df, lanes=(2, 3, 4, 5, 7, 8, 9)):
    out = []
    for k in lanes:
        sub = df[df.lane == k].sort_values("y_px")
        if len(sub) == 2:
            out.append(sub.volume.iloc[0] / sub.volume.iloc[1])
    return out


def test_touching_bands_split_by_peak_fit(tmp_path):
    """Resolved doublet 4:1 - the valley split gives the faint band the strong one's tail."""
    p = _doublet_gel(tmp_path, 14, 4)
    res = {}
    for dec in (False, True):
        ga = run(p, n_lanes=10, deconvolve=dec)
        df = ga.quantify()
        res[dec] = _ratios(df)
        assert len(res[dec]) == 7 and df["deconvolved"].any() == dec
    assert np.mean(res[True]) == pytest.approx(4.0, rel=0.03)
    assert abs(np.mean(res[False]) - 4) > 2 * abs(np.mean(res[True]) - 4)


@pytest.mark.parametrize("sep, ratio", [(10, 3), (9, 1.5), (14, 4)])
def test_split_shoulders(tmp_path, sep, ratio):
    p = _doublet_gel(tmp_path, sep, ratio)
    assert _ratios(run(p, n_lanes=10).quantify()) == [] or sep >= 14     # unresolved without the option
    ga = run(p, n_lanes=10, split_shoulders=True)
    df = ga.quantify()
    r = _ratios(df)
    assert len(r) == 7 and np.mean(r) == pytest.approx(ratio, rel=0.05)
    assert len(df[df.lane == 6]) == 1                                    # broad single band not split


def test_split_shoulders_skips_saturated(tmp_path):
    lanes = {k: [(500, 60000, 6.0)] for k in range(1, 9)}               # clipped flat tops
    ga = run(make_gel(str(tmp_path / "s.tif"), lanes=lanes, clip=30000), n_lanes=10, split_shoulders=True)
    assert all(sum(b.lane == k for b in ga.bands) == 1 for k in range(1, 9))


def test_noise_floor_keeps_faint_ladder_and_drops_noise(tmp_path):
    """Empty lanes of a noisy image produce noise 'bands'; the floor removes them, while a faint
    band far above the noise (but < 1 % of the bright samples) stays."""
    lad, _ = ladder(sizes=SIZES, amp=8000)
    lanes = {0: lad, 1: [(600, 60000 / 8, 5)], 2: [(600, 400, 4)]}         # faint but clean band in lane 3
    p = make_gel(str(tmp_path / "n.tif"), lanes=lanes, noise=60)
    on = run(p, n_lanes=10, noise_floor=True, min_snr=3)       # low SNR threshold: noise gets through
    off = run(p, n_lanes=10, noise_floor=False, min_snr=3)
    assert sum(b.lane == 2 for b in on.bands) == 1
    assert sum(b.lane == 0 for b in on.bands) == len(SIZES)
    assert sum(b.lane >= 3 for b in off.bands) >= 1 and sum(b.lane >= 3 for b in on.bands) == 0
    assert "noise floor" in on.log[-1]["result"]


def test_curved_lanes_traced_straight_lanes_left_alone(tmp_path):
    from .gallery import SCENARIOS, render
    by = {s.name: s for s in SCENARIOS}
    for name in ("agarose_clean", "agarose_low_signal", "agarose_uneven_light", "agarose_dirty"):
        p, _ = render(by[name], str(tmp_path))       # noise, vignetting, dust must not look curved
        ga = GelAnalysis(p, Params(n_lanes=12, curved_lanes=True))
        ga.detect_lanes()
        assert not any(l.bend for l in ga.lanes), name
    p, _ = render(by["agarose_curved_lanes"], str(tmp_path))
    ga = GelAnalysis(p, Params(n_lanes=12, curved_lanes=True))
    ga.detect_lanes()
    bent = [l for l in ga.lanes if l.bend]
    assert len(bent) >= 6
    right, left = ga.lanes[-1], ga.lanes[0]               # edge lanes drift ~25 px outward in the simulation
    assert right.bend and 12 < right.dx(right.y1) - right.dx(right.y0 + 60) < 40
    assert left.bend and -40 < left.dx(left.y1) - left.dx(left.y0 + 60) < -12
    # the profile follows the lane: a session round trip keeps the bend
    ga.detect_bands()
    s = tmp_path / "s.json"
    ga.save_session(str(s))
    gb = GelAnalysis.load_session(str(s))
    assert gb.lanes[-1].bend == ga.lanes[-1].bend
    assert np.allclose(gb.lane_profile(11)["corr"], ga.lane_profile(11)["corr"])
    # default (off): straight lanes
    gc = GelAnalysis(p, Params(n_lanes=12))
    gc.detect_lanes()
    assert not any(l.bend for l in gc.lanes)


@pytest.mark.parametrize("ext", [".tif", ".jpg"])
def test_cmyk_image_read_as_light_not_ink(tmp_path, ext):
    """Journal figure files are often CMYK; reading the ink channels as RGB inverted the image."""
    from PIL import Image

    from autogel.imageio import load_raw
    path = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), dtype=np.uint8, gradient=0, bg=20, noise=2)
    rgb = Image.open(path).convert("RGB")
    out = str(tmp_path / f"cmyk{ext}")
    rgb.convert("CMYK").save(out)
    assert Image.open(out).mode == "CMYK"
    a, _, _ = load_raw(out)
    ref = np.asarray(rgb.convert("L"), float)
    assert np.corrcoef(np.asarray(a, float).ravel(), ref.ravel())[0, 1] > 0.95


def test_very_large_image_is_reduced(tmp_path):
    """~36 MP image: analysed at 1/2 resolution (memory), relative volumes unchanged."""
    from autogel.imageio import MAX_MEGAPIXELS
    small = make_gel(str(tmp_path / "s.tif"), lanes=standard_lanes(), noise=10)
    a = tifffile.imread(small)
    big = str(tmp_path / "big.tif")
    tifffile.imwrite(big, np.kron(a, np.ones((5, 5), a.dtype)))      # 6000 x 5000 = 30 MP
    assert a.size * 25 / 1e6 > MAX_MEGAPIXELS
    ga = GelAnalysis(big, Params(n_lanes=10))
    assert ga.gel.reduced == 2 and ga.data.shape == (3000, 2500)
    assert any("Very large image" in w for w in ga.gel.warnings)
    ga.detect_lanes()
    ga.detect_bands()
    df = ga.quantify()
    ref = run(small, n_lanes=10).quantify()
    for k in (3, 5):
        r_big = df[df.lane == k].sort_values("y_px").volume.to_numpy()
        r_ref = ref[ref.lane == k].sort_values("y_px").volume.to_numpy()
        assert len(r_big) == len(r_ref) == 2
        assert r_big[0] / r_big[1] == pytest.approx(r_ref[0] / r_ref[1], rel=0.05)
    # an added image of the same (full) size is reduced the same way
    k = ga.add_image(big, name="copy")
    assert ga.gel.channels[k].shape == ga.data.shape


def _gel_in_frame(tmp_path):
    """Gel in the middle of a bright frame (tray / film edge) with lane-like bright stripes."""
    p = make_gel(str(tmp_path / "g.tif"), lanes=standard_lanes(), noise=20)
    g = tifffile.imread(p).astype(np.float32)
    H, W = g.shape
    big = np.full((H + 400, W + 600), 300, np.float32)       # dark tray around the gel
    big[200:200 + H, 300:300 + W] = g
    for x in (60, 120, 180, W + 360, W + 420, W + 480):      # stripes on the frame
        big[:, x:x + 30] += 6000
    out = str(tmp_path / "framed.tif")
    tifffile.imwrite(out, np.clip(big, 0, 65535).astype(np.uint16))
    return out, (200, 200 + H, 300, 300 + W)


def test_auto_crop_finds_the_gel_and_ignores_the_frame(tmp_path):
    from autogel.analysis import find_gel_region
    path, (y0, y1, x0, x1) = _gel_in_frame(tmp_path)
    ga = GelAnalysis(path, Params(n_lanes=10, auto_crop=True))
    ry0, ry1, rx0, rx1 = ga.region()
    assert abs(ry0 - y0) < 40 and abs(ry1 - y1) < 40 and abs(rx0 - x0) < 40 and abs(rx1 - x1) < 40
    ga.detect_lanes()
    assert all(x0 - 10 < l.x < x1 + 10 for l in ga.lanes)
    xs = np.array([l.x - x0 for l in ga.lanes])
    assert np.allclose(xs, 140 + 80 * np.arange(10), atol=6)
    # a gel that fills the image is not cropped
    full = make_gel(str(tmp_path / "full.tif"), lanes=standard_lanes())
    assert find_gel_region(tifffile.imread(full).astype(np.float32)) is None
    # a crop set by hand wins and is logged
    ga.set_crop([y0, y1, x0 + 80, x1])
    assert ga.region()[2] == x0 + 80 and ga.log[-1]["step"] == "set_crop"

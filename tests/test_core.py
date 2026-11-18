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

"""Every control in the napari panel must do something measurable.

Runs the real widget on a synthetic gel (known lanes/bands). File dialogs and message boxes
are answered automatically. Skipped if napari/Qt are not installed.
"""
import os
import time

import numpy as np
import pytest

napari = pytest.importorskip("napari")
pytest.importorskip("qtpy")

from napari.qt import get_qapp  # noqa: E402
from qtpy.QtWidgets import QFileDialog, QInputDialog, QMessageBox  # noqa: E402

from .synth import ladder, make_gel, y_for_size  # noqa: E402

SIZES = [10000, 6000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 100]


def pump(t=1.5):
    app = get_qapp()
    end = time.time() + t
    while time.time() < end:
        app.processEvents()
        time.sleep(0.01)


@pytest.fixture
def ui(tmp_path, monkeypatch):
    monkeypatch.setenv("AUTOGEL_HOME", str(tmp_path / "home"))
    from qtpy.QtCore import QSettings
    QSettings("autogel", "autogel").setValue("mode", "easy")
    msgs = []
    for name in ("information", "warning", "critical"):
        monkeypatch.setattr(QMessageBox, name, lambda *a, _n=name, **k: msgs.append((_n, a[2] if len(a) > 2 else "")))
    get_qapp()
    from autogel.napari_widget import AutogelWidget
    lad, _ = ladder(sizes=SIZES)
    lanes = {0: lad, 9: lad}
    for k in range(1, 9):
        lanes[k] = [(y_for_size(2500), 6000, 5), (y_for_size(800), 3000, 5)]
    lanes[5] = lanes[5] + [(y_for_size(5000), 150, 5)]     # faint band (S/N ~ 16)
    path = make_gel(str(tmp_path / "gel.tif"), lanes=lanes)
    v = napari.Viewer(show=False)
    w = AutogelWidget(v)
    w.open_image(path)
    yield w, v, msgs, tmp_path
    w.shutdown()
    pump(0.2)
    v.close()
    pump(0.2)


def ready(w, n=10):
    w.sp_lanes.setValue(n)
    w.detect_lanes()
    w.detect_bands()
    w.ed_marker_lanes.setText("1, 10")
    w.cb_preset.setCurrentText("custom")
    w.ed_sizes.setText(", ".join(map(str, SIZES)))
    w.ed_masses.setText(", ".join(["10"] * len(SIZES)))
    w.quantify()


# ----------------------------------------------------------------- panel
def test_mode_switch_hides_expert_controls_and_persists(ui):
    w, v, msgs, _ = ui
    from qtpy.QtCore import QSettings
    assert not w.cb_bg.isVisibleTo(w) and w.cb_preset.isVisibleTo(w)
    w.btn_expert.setChecked(True)
    assert w.cb_bg.isVisibleTo(w) and w.ed_sizes.isVisibleTo(w)
    assert QSettings("autogel", "autogel").value("mode") == "expert"
    w.btn_easy.setChecked(True)
    assert not w.sp_snr.isVisibleTo(w)


def test_status_line_and_log_panel(ui):
    w, v, msgs, _ = ui
    assert "Loaded gel.tif" in w.lbl_status.text()
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    assert "Detected 10 lanes" in w.lbl_status.text()
    assert not w.txt_log.isVisibleTo(w)
    w.btn_log.setChecked(True)
    assert w.txt_log.isVisibleTo(w)
    assert "Detected 10 lanes" in w.txt_log.toPlainText()


# ------------------------------------------------------------ orientation
def test_rotate_and_flips(ui):
    w, v, msgs, _ = ui
    a0 = w.ga.data.copy()
    w._flip_lr()
    assert np.allclose(w.ga.data, a0[:, ::-1])
    w._flip_lr()
    w._flip_ud()
    assert np.allclose(w.ga.data, a0[::-1, :])
    w._flip_ud()
    w._rot90()
    assert w.ga.data.shape == a0.shape[::-1] and w.ga.gel.orientation["rot90_k"] % 4 != 0
    assert any(e["step"] == "orientation_change" for e in w.ga.log)


def test_invert_and_skew(ui):
    w, v, msgs, _ = ui
    a0 = w.ga.data.copy()
    w.chk_invert.setChecked(True)
    assert w.ga.gel.orientation["invert"] is True
    assert np.corrcoef(w.ga.data.ravel()[::97], a0.ravel()[::97])[0, 1] < -0.9
    w.chk_invert.setChecked(False)
    w.sp_angle.setValue(1.5)
    w.sp_angle.editingFinished.emit()
    assert w.ga.gel.orientation["angle"] == pytest.approx(1.5)


# ------------------------------------------------------------------ lanes
def test_lane_count_given_and_asked(ui, monkeypatch):
    w, v, msgs, _ = ui
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    assert len(w.ga.lanes) == 10 and len(v.layers["lanes"].data) == 10
    w.sp_lanes.setValue(0)       # 'ask': user corrects the estimate to 9
    monkeypatch.setattr(QInputDialog, "getInt", lambda *a, **k: (9, True))
    w.detect_lanes()
    assert len(w.ga.lanes) == 9 and w.sp_lanes.value() == 9


def test_lane_names(ui):
    w, v, msgs, _ = ui
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    w.ed_names.setText("M, ctrl, s1")
    w.ed_names.editingFinished.emit()
    assert [l.name for l in w.ga.lanes[:4]] == ["M", "ctrl", "s1", "4"]
    assert list(v.layers["lanes"].features["name"][:2]) == ["M", "ctrl"]


# ------------------------------------------------------------------ bands
def test_band_detection_parameters_have_effect(ui):
    w, v, msgs, _ = ui
    ready(w)
    n0 = len(w.ga.bands)
    vol0 = w.df.volume.sum()
    w.sp_snr.setValue(60)
    w.detect_bands()
    assert len(w.ga.bands) == n0 - 1               # stricter S/N drops the faint band
    w.sp_snr.setValue(6)
    w.cb_bg.setCurrentText("none")
    w.detect_bands()
    assert w.df.volume.sum() > vol0 * 1.2           # no background subtraction -> larger volumes
    w.cb_bg.setCurrentText("lane_edges")
    w.detect_bands()
    assert w.ga.params.bg_method == "lane_edges" and len(w.ga.bands) > 0
    w.cb_bg.setCurrentText("rolling_ball")
    w.sp_radius.setValue(10)
    w.detect_bands()
    assert w.ga.params.bg_radius == 10
    w.sp_gelrel.setValue(0.1)                      # gel-relative floor drops the faint band
    w.sp_snr.setValue(6)
    w.cb_bg.setCurrentText("rolling_ball")
    w.detect_bands()
    assert w.ga.params.min_gel_rel == pytest.approx(0.1)
    assert any("below 10%" in e["result"] for e in w.ga.log if e["step"] == "detect_bands")
    w.sp_gelrel.setValue(0)
    w.sp_cov.setValue(0.99)
    w.detect_bands()
    n_strict = len(w.ga.bands)
    w.sp_cov.setValue(0.0)
    w.detect_bands()
    assert len(w.ga.bands) >= n_strict


def test_redetect_selected_lane(ui):
    w, v, msgs, _ = ui
    ready(w)
    keep = [b for b in w.ga.bands if b.lane != 3]
    w.ga.set_bands(keep)
    w._refresh_layers()
    v.layers["lanes"].selected_data = {3}
    w._redetect_selected()
    assert sum(b.lane == 3 for b in w.ga.bands) == 2
    v.layers["lanes"].selected_data = set()
    w._redetect_selected()
    assert any("Select a lane" in m[1] for m in msgs)


def test_viewer_edits_delete_add_move(ui):
    w, v, msgs, _ = ui
    ready(w)
    from autogel.napari_widget import _rect
    n = len(w.ga.bands)
    bl = v.layers["bands"]
    bl.selected_data = {0}
    bl.remove_selected()
    pump()
    assert len(w.ga.bands) == n - 1
    lane = w.ga.lanes[4]
    bl.add_rectangles([_rect(1000, lane.x0 + 2, 1015, lane.x1 - 2)])
    pump()
    assert len(w.ga.bands) == n and any(b.source == "user" for b in w.ga.bands)
    ll = v.layers["lanes"]
    d = list(ll.data)
    d[2] = d[2] + np.array([0, 6])
    ll.data = d
    pump()
    assert w.ga.lanes[2].x == pytest.approx(140 + 2 * 80 + 6, abs=4)


# ------------------------------------------------------------ calibration
def test_ladder_preset_fills_fields_and_info(ui):
    w, v, msgs, _ = ui
    from autogel.markers import DEFAULT_PRESET
    w.cb_preset.setCurrentText("GeneRuler 1 kb (SM0311)")
    assert w.ed_sizes.text().startswith("10000, 8000, 6000")
    assert w.sp_ref.value() == pytest.approx(0.5)
    assert "checked against vendor data" in w.lbl_ladder.text()
    w.cb_preset.setCurrentText(DEFAULT_PRESET)
    assert w.sp_ref.value() == pytest.approx(1.0) and "10002" in w.ed_sizes.text()


def test_calibration_sizes_masses_loads(ui):
    w, v, msgs, _ = ui
    ready(w)
    assert w.table.rowCount() == len(w.df) > 0
    s = w.df[w.df.lane == 3].sort_values("y_px")["size_bp"].tolist()
    assert s[0] == pytest.approx(2500, rel=0.05)
    m1 = w.df[w.df.lane == 3]["est_mass_ng"].iloc[0]
    w.sp_ref.setValue(1.0)
    w.sp_load.setValue(2.0)        # twice the reference amount loaded -> masses double
    w.quantify()
    assert w.df[w.df.lane == 3]["est_mass_ng"].iloc[0] == pytest.approx(2 * m1, rel=2e-3)
    w.ed_loads.setText("3:5, nonsense, 4:x")
    w.quantify()
    assert w.df[w.df.lane == 3]["conc_ng_per_ul"].iloc[0] == pytest.approx(2 * m1 / 5, rel=1e-3)
    assert "conc_ng_per_ul" in w.df and np.isnan(w.df[w.df.lane == 4]["conc_ng_per_ul"]).all()


def test_save_own_ladder(ui, monkeypatch):
    w, v, msgs, tmp = ui
    w.cb_preset.setCurrentText("custom")
    w.ed_sizes.setText("3000, 1000, 500")
    w.ed_masses.setText("30, 20, 10")
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("lab ladder", True))
    monkeypatch.setattr(QInputDialog, "getItem", lambda *a, **k: ("bp", True))
    w._save_ladder()
    items = [w.cb_preset.itemText(i) for i in range(w.cb_preset.count())]
    assert "lab ladder" in items and w.cb_preset.currentText() == "lab ladder"
    assert os.path.exists(tmp / "home" / "ladders.json")
    from autogel import markers
    markers.PRESETS.pop("lab ladder", None)


def test_bad_marker_input_reports_instead_of_crashing(ui):
    w, v, msgs, _ = ui
    ready(w)
    w.ed_sizes.setText("10000, abc")
    w.quantify()
    assert any(m[0] == "critical" and "Could not read" in m[1] for m in msgs)


# ---------------------------------------------------------------- display
def test_palette_and_label_all(ui):
    w, v, msgs, _ = ui
    ready(w)
    edge0 = np.array(v.layers["lanes"].edge_color[1])
    w.cb_palette.setCurrentText("classic")
    assert not np.allclose(np.array(v.layers["lanes"].edge_color[1]), edge0)
    n_marker_only = len(v.layers["sizes left"].data) + len(v.layers["sizes right"].data)
    v.scene.camera.zoom = 5            # zoom in: no decluttering
    w.chk_all_labels.setChecked(True)
    pump(0.3)
    n_all = len(v.layers["sizes left"].data) + len(v.layers["sizes right"].data)
    assert n_all > n_marker_only


def test_profile_navigation_and_hover(ui):
    w, v, msgs, _ = ui
    ready(w)
    p = w.profile
    p.step(+1)
    assert p.lane == 1 and w._focus_lane == 1
    p.sp_lane.setValue(5)
    assert p.lane == 4
    for _ in range(20):
        p.step(+1)
    assert p.lane == 9 and not p.btn_next.isEnabled()
    p._on_key(type("K", (), {"key": "left"})())
    assert p.lane == 8
    k = p._band_ids[0]
    b = w.ga.bands[k]
    p._on_move(type("E", (), {"inaxes": p.ax, "xdata": (b.y0 + b.y1) / 2})())
    assert p._hover == k and len(v.layers["highlight"].data) == 2


# ----------------------------------------------------------------- output
def test_export_and_session_buttons(ui, monkeypatch):
    w, v, msgs, tmp = ui
    ready(w)
    out = tmp / "out"
    out.mkdir()
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: str(out))
    w.export()
    names = os.listdir(out)
    for suffix in ("_report.md", "_log.jsonl", "_bands.csv", "_annotated.png", "_session.json"):
        assert "gel" + suffix in names
    assert "Exported" in w.lbl_status.text()
    sess = str(tmp / "s.json")
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (sess, ""))
    w.save_session()
    assert os.path.exists(sess)
    n = len(w.ga.bands)
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (sess, ""))
    w._open_session_dialog()
    assert len(w.ga.bands) == n and "Session restored" in w.lbl_status.text()
    assert w.ed_marker_lanes.text() == "1, 10"


def test_size_range_control(ui):
    w, v, msgs, _ = ui
    ready(w)
    n_cal = len(w.ga.calibrations[0].ys)
    w.cb_range.setCurrentText("agarose 1.5 %")            # 200-3000 bp
    assert (w.sp_rmin.value(), w.sp_rmax.value()) == (200, 3000)
    assert len(w.ga.calibrations[0].ys) < n_cal
    assert w.df[w.df.lane == 3]["size_in_range"].all()    # 2500 and 800 bp are inside
    w.sp_rmax.setValue(2000)                              # typing a value -> custom
    assert w.cb_range.currentText() == "custom"
    assert not w.df[w.df.lane == 3].sort_values("y_px")["size_in_range"].iloc[0]  # 2500 now outside
    assert w.ga.params.size_max == 2000
    w.cb_range.setCurrentText("automatic (whole ladder)")
    assert w.ga.params.size_min is None and len(w.ga.calibrations[0].ys) == n_cal
    assert any(e["step"] == "set_size_range" for e in w.ga.log)


def test_two_gel_rows_and_split_line(ui):
    w, v, msgs, tmp = ui
    from .gallery import SCENARIOS, render
    sc = {s.name: s for s in SCENARIOS}["agarose_two_combs"]
    path, _ = render(sc, str(tmp))
    w.open_image(path)
    w.sp_tiers.setValue(2)
    w.sp_lanes.setValue(12)
    w.detect_lanes()
    assert len(w.ga.lanes) == 24 and {l.tier for l in w.ga.lanes} == {0, 1}
    assert w.ga.lanes[12].name == "B1"
    lines = v.layers["gel rows"].data
    assert len(lines) == 1 and 880 < np.mean(lines[0][:, 0]) < 990
    # drag the split line down by 20 px -> lanes re-detected for the new rows
    v.layers["gel rows"].data = [lines[0] + np.array([20, 0])]
    pump()
    assert any(e["step"] == "set_tiers" and "moved" in str(e["params"]) for e in w.ga.log)
    assert len(w.ga.lanes) == 24


def test_batch_button(ui, tmp_path):
    w, v, msgs, tmp = ui
    ready(w)
    from .synth import make_gel
    lad, _ = ladder(sizes=SIZES)
    imgs = [make_gel(str(tmp_path / f"b{i}.tif"), lanes={0: lad, 9: lad, 4: [(y_for_size(2500), 6000, 5)]},
                     seed=i) for i in range(2)]
    bad = tmp_path / "bad.tif"
    bad.write_bytes(b"nope")
    out = tmp_path / "batch"
    summary = w.batch(images=[*imgs, str(bad)], out_dir=str(out))
    assert list(summary.status) == ["ok", "ok", "error"]
    assert (summary.lanes.iloc[:2] == 10).all()
    for f in ("batch_report.md", "batch_bands.csv", "batch_summary.csv", "batch_autogel.xlsx"):
        assert (out / f).exists()
    import pandas as pd
    bands = pd.read_csv(out / "batch_bands.csv")
    s = bands[(bands.gel == "b0.tif") & (bands.lane == 5)]["size_bp"].iloc[0]
    assert s == pytest.approx(2500, rel=0.03)
    assert any("Batch finished" in m[1] for m in msgs)

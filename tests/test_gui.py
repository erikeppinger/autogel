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
    # a wrong automatic orientation must be fixable in Easy mode too
    from qtpy.QtWidgets import QPushButton
    btns = {b.text(): b for b in w.findChildren(QPushButton)}
    assert all(btns[t].isVisibleTo(w) for t in ("Rotate 90°", "Flip ↔", "Flip ↕"))
    assert not w.sp_angle.isVisibleTo(w) and not w.chk_invert.isVisibleTo(w)


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


def test_flip_reverses_skew(ui):
    """A mirrored image tilts the other way: Flip must negate the skew, not keep it (a 2.5°
    tilted gel was left 5° crooked before 0.3.1)."""
    w, v, msgs, tmp = ui
    from .gallery import SCENARIOS, render
    path, _ = render({s.name: s for s in SCENARIOS}["agarose_tilted"], str(tmp))
    w.open_image(path)
    ang = w.ga.gel.orientation["angle"]
    assert abs(ang) > 1
    w._flip_lr()
    assert w.ga.gel.orientation["angle"] == pytest.approx(-ang)
    assert w.sp_angle.value() == pytest.approx(-ang)
    w._flip_ud()
    assert w.ga.gel.orientation["angle"] == pytest.approx(ang)
    w._rot90()                                  # rotation keeps the tilt, re-measured
    assert w.ga.gel.orientation["angle"] == pytest.approx(ang, abs=0.3)


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


def test_ladder_manager_dialog(ui, monkeypatch, tmp_path):
    w, v, msgs, tmp = ui
    from autogel import markers
    from autogel.ladder_dialog import LadderDialog
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    d = LadderDialog(w, select=markers.DEFAULT_PRESET)
    assert d.ed_name.text() == markers.DEFAULT_PRESET and "built-in" in d.lbl_origin.text()
    # editing a built-in ladder saves a copy, never the original
    d.ed_note.setText("our lot")
    d._save()
    assert markers.DEFAULT_PRESET + " (my copy)" in markers.PRESETS
    assert markers.PRESETS[markers.DEFAULT_PRESET]["verified"]
    # new ladder; sizes in any order, masses stay with them
    d._new()
    d.ed_name.setText("lab ladder")
    d.ed_sizes.setText("500 3000 1000")
    d.ed_masses.setText("10, 30, 20")
    d.sp_ref.setValue(0.06)
    assert d._save() and markers.PRESETS["lab ladder"]["masses"] == [30, 20, 10]
    assert markers.PRESETS["lab ladder"]["reference_load"]["amount"] == pytest.approx(0.06)
    d.ed_sizes.setText("1000, abc")
    assert not d._save() and "Not saved" in d.lbl_msg.text()
    # export, delete, import back
    out = tmp_path / "lab.json"
    d._show("lab ladder")
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    d._export()
    assert out.exists()
    d._delete()
    assert "lab ladder" not in markers.PRESETS
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(out), ""))
    d._import()
    assert "lab ladder" in markers.PRESETS and "1 added" in d.lbl_msg.text()
    d._delete()                                     # deleting a built-in is refused
    d._show(markers.DEFAULT_PRESET)
    d._delete()
    assert markers.DEFAULT_PRESET in markers.PRESETS and "own ladders" in d.lbl_msg.text()
    # sharing opens a pre-filled issue (no network access here: the URL is captured)
    opened = []
    from qtpy.QtGui import QDesktopServices
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda u: opened.append(u.toString()))
    d.ed_source.setText("")
    d._share()
    assert not opened and "datasheet" in d.lbl_msg.text()
    d.ed_source.setText("https://example.org/ds")
    d._share()
    assert opened and "issues/new" in opened[0] and "template=ladder.yml" in opened[0]
    # update failure keeps the old library
    monkeypatch.setattr(markers, "update_library", lambda *a, **k: (_ for _ in ()).throw(OSError("offline")))
    d._update()
    assert "old library kept" in d.lbl_msg.text() and markers.DEFAULT_PRESET in markers.PRESETS
    # "use" with unsaved edits saves them and selects the ladder in the panel
    d._new()
    d.ed_name.setText("quick")
    d.ed_sizes.setText("2000, 1000, 500")
    d._use()
    assert d.chosen == "quick"
    w._refresh_presets(d.chosen)
    assert w.cb_preset.currentText() == "quick" and w.ed_sizes.text() == "2000, 1000, 500"
    for n in ("quick", "lab ladder", markers.DEFAULT_PRESET + " (my copy)"):
        markers.delete_user_ladder(n)


def test_suggest_ladder_button(ui, monkeypatch):
    w, v, msgs, _ = ui
    ready(w)
    w.ed_marker_lanes.setText("1")
    seen = {}

    def pick(parent, title, label, items, *a, **k):
        seen["items"] = items
        return items[0], True
    monkeypatch.setattr(QInputDialog, "getItem", pick)
    w._suggest_ladder()
    assert seen["items"] and any(e["step"] == "suggest_ladder" for e in w.ga.log)
    assert w.cb_preset.currentText() in [i.split("   (")[0] for i in seen["items"]]


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
    lines = v.layers["gel splits"].data
    assert len(lines) == 1 and 880 < np.mean(lines[0][:, 0]) < 990
    # drag the split line down by 20 px -> lanes re-detected for the new rows
    v.layers["gel splits"].data = [lines[0] + np.array([20, 0])]
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
    # Cancel in the progress window (as pressed by hand in the GUI test on the lab gels)
    from qtpy.QtCore import QTimer
    from qtpy.QtWidgets import QApplication, QProgressDialog

    def press_cancel():
        dlg = next((d for d in QApplication.topLevelWidgets()
                    if isinstance(d, QProgressDialog) and d.isVisible()), None)
        if dlg is None:
            QTimer.singleShot(50, press_cancel)
        else:
            dlg.cancel()
    QTimer.singleShot(200, press_cancel)
    summary = w.batch(images=imgs * 3, out_dir=str(tmp_path / "cancelled"))
    assert summary.attrs["cancelled"] and len(summary) < 6
    assert any(f"Batch cancelled after {len(summary)} of 6 gel(s)" in m[1] for m in msgs)
    assert w.lbl_status.text().startswith("✔ Batch cancelled")
    # the real button: Qt passes "checked" to the slot, which once arrived as images=False
    # and made the button silently do nothing
    asked = []
    monkey = pytest.MonkeyPatch()
    monkey.setattr(QFileDialog, "getOpenFileNames",
                   staticmethod(lambda *a, **k: (asked.append("files"), (imgs[:1], ""))[1]))
    monkey.setattr(QFileDialog, "getExistingDirectory",
                   staticmethod(lambda *a, **k: (asked.append("folder"), str(tmp_path / "button"))[1]))
    try:
        w.btn_batch.click()
    finally:
        monkey.undo()
    assert asked == ["files", "folder"]
    assert (tmp_path / "button" / "batch_report.md").exists()


def test_sample_sheet_buttons(ui, tmp_path):
    w, v, msgs, tmp = ui
    ready(w)
    import pandas as pd

    from .synth import make_gel
    lad, _ = ladder(sizes=SIZES)
    imgs = [make_gel(str(tmp_path / "a.tif"), lanes={0: lad, 9: lad}),
            make_gel(str(tmp_path / "b.tif"), n_lanes=6, lanes={0: lad, 3: [(y_for_size(2500), 6000, 5)]}, seed=1)]
    sheet = str(tmp_path / "sheet.xlsx")
    monkey = pytest.MonkeyPatch()
    monkey.setattr(QFileDialog, "getOpenFileNames", staticmethod(lambda *a, **k: (imgs, "")))
    monkey.setattr(QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (sheet[:-5], "")))
    try:
        w.btn_make_sheet.click()                        # ".xlsx" is added when missing
    finally:
        monkey.undo()
    df = pd.read_excel(sheet, dtype=object)
    assert list(df.image) == ["a.tif", "b.tif"] and df.marker_lanes[0] == "1, 10"
    assert any("Sample sheet saved" in m[1] for m in msgs)
    df.loc[1, ["lanes", "marker_lanes", "lane_names"]] = [6, "1", "L, x, y, s, z, e"]
    df.to_excel(sheet, index=False)
    monkey = pytest.MonkeyPatch()
    monkey.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (sheet, "")))
    monkey.setattr(QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(tmp_path / "out")))
    try:
        w.btn_sheet_batch.click()
    finally:
        monkey.undo()
    s = pd.read_csv(tmp_path / "out" / "batch_summary.csv")
    assert s.lanes.tolist() == [10, 6]
    assert any("Batch finished: 2 gel(s)" in m[1] for m in msgs)
    # a sheet with a mistake: a plain message, no batch
    pd.DataFrame({"image": ["a.tif"], "ladder": ["no such ladder"]}).to_excel(sheet, index=False)
    w.batch_from_sheet(sheet, str(tmp_path / "out2"))
    assert msgs[-1][0] == "critical" and "row 2: ladder 'no such ladder'" in msgs[-1][1]
    assert not (tmp_path / "out2").exists()


def test_next_band_to_check(ui):
    w, v, msgs, tmp = ui
    ready(w)
    assert w.btn_next_check.isVisibleTo(w)            # in Easy mode too
    w.df = w.df.assign(band_probability=1.0)
    w.btn_next_check.click()
    assert "No bands left to check" in w.lbl_status.text()
    # two unsure bands: the more doubtful one comes first, then they cycle
    samples = [k for k, b in enumerate(w.ga.bands) if b.lane not in (0, 9)]
    i, j = samples[2], samples[6]
    k_marker = next(k for k, b in enumerate(w.ga.bands) if b.lane == 0)
    w.df.loc[k_marker, "band_probability"] = 0.001     # ladder lanes: checked by the calibration
    w.df.loc[i, "band_probability"] = 0.05
    w.df.loc[j, "band_probability"] = 0.01
    w._refresh_layers()
    assert w._doubtful() == [j, i]
    bands = v.layers["bands"]
    for expect, nth in ((j, 1), (i, 2), (j, 1)):
        w.btn_next_check.click()
        assert bands.selected_data == {expect}
        assert w._focus_lane == w.ga.bands[expect].lane
        assert f"Band {nth} of 2 to check" in w.lbl_status.text()
    # grey colour follows the same rule
    from matplotlib.colors import to_rgba

    from autogel.napari_widget import DOUBT_BAND
    assert np.allclose(bands.edge_color[j], to_rgba(DOUBT_BAND))
    assert not np.allclose(bands.edge_color[0], to_rgba(DOUBT_BAND))
    # delete the band just shown (as with the Delete key): the next click must not skip one
    w.btn_next_check.click()
    assert bands.selected_data == {i}
    w.btn_next_check.click()
    assert bands.selected_data == {j}
    n = len(w.ga.bands)
    bands.remove_selected()
    pump(1.5)                                          # edits apply after a short pause
    assert len(w.ga.bands) == n - 1
    w.df = w.df.assign(band_probability=1.0)           # the edit recomputed the table
    w.df.loc[i, "band_probability"] = 0.05             # i < j: its index did not shift
    assert i < j
    w._refresh_layers()
    w.btn_next_check.click()
    assert bands.selected_data == {i}
    assert "Band 1 of 1 to check" in w.lbl_status.text()


def test_gels_side_by_side_and_split_line(ui):
    w, v, msgs, tmp = ui
    from .gallery import SCENARIOS, render
    sc = {s.name: s for s in SCENARIOS}["agarose_side_by_side"]
    path, _ = render(sc, str(tmp))
    w.open_image(path)
    w.sp_cols.setValue(2)
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    assert len(w.ga.lanes) == 20 and {l.column for l in w.ga.lanes} == {0, 1}
    assert [w.ga.lanes[0].name, w.ga.lanes[10].name] == ["A1", "B1"]
    lines = v.layers["gel splits"].data
    assert len(lines) == 1 and 880 < np.mean(lines[0][:, 1]) < 980
    v.layers["gel splits"].data = [lines[0] + np.array([0, 10])]     # drag the split right
    pump()
    assert any(e["step"] == "set_columns" and "moved" in str(e["params"]) for e in w.ga.log)
    assert len(w.ga.lanes) == 20
    # markers of gel A must not calibrate gel B: give markers only in gel A
    w.ed_marker_lanes.setText("1, 10")
    w.cb_preset.setCurrentText("custom")
    from .gallery import DNA_LADDER
    w.ed_sizes.setText(", ".join(map(str, DNA_LADDER)))
    w.quantify()
    assert w.df[w.df.lane > 10]["size_bp"].isna().all()


def test_curved_lanes_drawn_and_editable(ui):
    w, v, msgs, tmp = ui
    from .gallery import SCENARIOS, render
    path, _ = render({s.name: s for s in SCENARIOS}["agarose_curved_lanes"], str(tmp))
    w.open_image(path)
    w.chk_curved.setChecked(True)
    w.sp_lanes.setValue(12)
    w.detect_lanes()
    w.detect_bands()
    layer = v.layers["lanes"]
    kinds = list(layer.shape_type)
    assert "polygon" in kinds and kinds.count("polygon") == sum(bool(l.bend) for l in w.ga.lanes)
    n_log = len(w.ga.log)
    layer.data = list(layer.data)            # data event without a real change
    pump()
    assert not any(e["step"] == "set_lanes" for e in w.ga.log[n_log:])
    # bands sit on the bent lane
    k = max(range(len(w.ga.lanes)), key=lambda i: abs(w.ga.lanes[i].dx(w.ga.lanes[i].y1)))
    lane = w.ga.lanes[k]
    b = max((b for b in w.ga.bands if b.lane == k), key=lambda b: b.y)
    box = [np.asarray(s) for s in v.layers["bands"].data][w.ga.bands.index(b)]
    assert abs(box[:, 1].mean() - lane.x_at(b.y)) < 1.5
    # drag the curved lane 6 px right: still curved, same name, edit logged
    data = [np.asarray(s, float) for s in layer.data]
    data[k] = data[k] + np.array([0, 6])
    layer.data = data
    pump()
    new = w.ga.lanes[k]
    assert any(e["step"] == "set_lanes" for e in w.ga.log[n_log:])
    assert new.bend and new.name == lane.name
    assert new.x_at(new.y1) == pytest.approx(lane.x_at(lane.y1) + 6, abs=1.0)


def test_new_band_and_mass_options_reach_the_analysis(ui, tmp_path, monkeypatch):
    w, v, msgs, _ = ui
    w.btn_expert.setChecked(True)
    for chk in (w.chk_floor, w.chk_deconv, w.chk_shoulder, w.chk_curved):
        assert chk.isVisibleTo(w)
    w.chk_floor.setChecked(False)
    w.chk_shoulder.setChecked(True)
    w.chk_curved.setChecked(False)          # default off; checked elsewhere
    ready(w)
    p = w.ga.params
    assert (p.noise_floor, p.split_shoulders, p.curved_lanes) == (False, True, False)
    w.chk_deconv.setChecked(False)
    w.cb_mass.setCurrentText("one slope through origin")
    w.quantify()
    assert p.deconvolve is False and p.mass_model == "linear"
    assert w.ga.mass_fit["model"] == "linear" and "through origin" in w.lbl_cal.text()
    w.cb_mass.setCurrentText("size-dependent")
    w.quantify()
    assert w.ga.mass_fit["model"] in ("size", "linear") and "typical error" in w.lbl_cal.text()
    # a saved session restores the switches
    s = tmp_path / "s.json"
    w.ga.save_session(str(s))
    w.chk_shoulder.setChecked(False)
    w.chk_floor.setChecked(True)
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(s), ""))
    w._open_session_dialog()
    assert w.chk_shoulder.isChecked() and not w.chk_floor.isChecked() and not w.chk_curved.isChecked()
    assert w.cb_mass.currentText() == "size-dependent"


def test_western_blot_workflow(ui, monkeypatch, tmp_path):
    """Chemiluminescence image + white-light marker image, loading-control normalisation."""
    w, v, msgs, _ = ui
    from .test_compare import AMOUNT, _western
    chemi, marker = _western(tmp_path)
    w.open_image(chemi)
    assert w.cb_channel.count() == 1 and not w.cb_marker_ch.isEnabled()      # plain grey image
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    w.detect_bands()
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (marker, ""))
    w._add_image_dialog()
    assert [w.cb_marker_ch.itemText(i) for i in range(w.cb_marker_ch.count())] == \
        ["same image", "1: image", "2: marker"]
    assert w.cb_channel.itemText(w.cb_channel.count() - 1) == "all channels"
    w.cb_preset.setCurrentText("Bio-Rad Precision Plus Protein Dual Color (1610374)")
    w.ed_marker_lanes.setText("1")
    w.cb_marker_ch.setCurrentIndex(2)
    assert w.ga.params.marker_channel == 1 and sum(b.lane == 0 for b in w.ga.bands) == 10
    w.cb_norm.setCurrentText("loading-control band")
    w.sp_lc.setValue(42)
    w.sp_ref_lane.setValue(2)
    w.quantify()
    tgt = w.df[(~w.df.is_marker) & w.df.size_kDa.between(50, 60)]
    fc = dict(zip(tgt.lane, tgt.fold_change))
    for k, amount in AMOUNT.items():
        assert fc[k + 1] == pytest.approx(amount, rel=0.08), k
    # switching the analysed channel shows that image and keeps lanes/bands
    n_bands = len(w.ga.bands)
    w.cb_channel.setCurrentIndex(2)
    assert w.ga.params.channel == 1 and np.allclose(v.layers["gel"].data, w.ga.gel.channels[1])
    assert len(w.ga.bands) == n_bands
    w.cb_channel.setCurrentIndex(0)
    # the session restores images and choices
    s = tmp_path / "wb.json"
    w.ga.save_session(str(s))
    w.cb_norm.setCurrentText("none")
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(s), ""))
    w._open_session_dialog()
    assert w.cb_norm.currentText() == "loading-control band" and w.sp_lc.value() == 42
    assert w.cb_marker_ch.currentData() == 1 and w.sp_ref_lane.value() == 2


def test_rna_checkbox_shows_ratio(ui, tmp_path):
    w, v, msgs, _ = ui
    from .synth import make_gel
    from .test_compare import RNA_LAD, y_nt
    lanes = {0: [(y_nt(s), 6000, 4) for s in RNA_LAD]}
    lanes.update({k: [(y_nt(5000), 8000, 5), (y_nt(1900), 4200, 5)] for k in range(1, 9)})
    w.open_image(make_gel(str(tmp_path / "rna.tif"), lanes=lanes, noise=30))
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    w.detect_bands()
    w.cb_preset.setCurrentText("NEB ssRNA Ladder (N0362)")
    w.ed_marker_lanes.setText("1")
    w.chk_rna.setChecked(True)
    w.quantify()
    assert "RNA 28S/18S" in w.lbl_cal.text()
    assert "rna_28s_18s" in w.ga.lanes_table()


def test_auto_crop_checkbox_and_region_drag(ui, tmp_path):
    w, v, msgs, _ = ui
    from .test_core import _gel_in_frame
    path, (y0, y1, x0, x1) = _gel_in_frame(tmp_path)
    w.open_image(path)
    w.btn_expert.setChecked(True)
    w.chk_crop.setChecked(True)
    w.sp_lanes.setValue(10)
    w.detect_lanes()
    assert w.ga.params.auto_crop and len(v.layers["analysed region"].data) == 1
    assert all(x0 - 10 < l.x < x1 + 10 for l in w.ga.lanes)
    rect = np.asarray(v.layers["analysed region"].data[0], float)
    v.layers["analysed region"].data = [rect + np.array([0, 40])]      # drag 40 px to the right
    pump()
    assert w.ga.params.crop is not None and abs(w.ga.region()[2] - (rect[:, 1].min() + 40)) < 2
    assert any(e["step"] == "set_crop" for e in w.ga.log)


def test_ai_button_and_checker_spin(ui, monkeypatch):
    w, v, msgs, _ = ui
    from autogel import ai
    ready(w)
    monkeypatch.setattr(ai, "available", lambda: (False, "the AI extra is not installed (pip install autogel[ai])"))
    w.ai_second_opinion()
    assert any("pip install autogel[ai]" in m[1] for m in msgs)
    from .test_models import _fake_network
    _fake_network(monkeypatch, w.ga, keep_lanes=set(range(10)))
    n0 = len(w.ga.bands)
    w.ai_second_opinion()
    assert len(w.ga.bands) >= n0 and "ai_support" in w.df and len(w.df) == len(w.ga.bands)
    assert len(v.layers["bands"].data) == len(w.ga.bands)
    assert any(e["step"] == "ai_review" for e in w.ga.log)
    w.btn_expert.setChecked(True)
    assert w.sp_checker.isVisibleTo(w)
    w.sp_checker.setValue(0.3)
    w.detect_bands()
    assert w.ga.params.checker_filter == pytest.approx(0.3)


def test_help_menu_guide_and_tour(ui, monkeypatch):
    w, v, msgs, _ = ui
    from autogel import help_dialog
    opened = []
    monkeypatch.setattr(help_dialog.QDesktopServices, "openUrl", lambda url: opened.append(url.toString()) or True)
    acts = {a.text(): a for a in w.btn_help.menu().actions() if a.text()}
    assert {"Quick tour", "User guide  (F1)", "Full manual (online)", "Report a problem or idea (online)"} <= set(acts)
    acts["Full manual (online)"].trigger()
    acts["Report a problem or idea (online)"].trigger()
    assert any("USER_MANUAL" in u for u in opened) and any("issues" in u for u in opened)
    # guide: offline, sections, search
    w.show_guide()
    g = w._guide
    assert g.isVisible() and g.toc.count() >= 8
    g.search.setText("saturated")
    assert 0 < g.toc.count() < len(g.sections)
    g.search.setText("no such word xyz")
    assert g.toc.count() == 0
    g.search.clear()
    w.show_guide("Step 4: Marker (ladder)")
    assert g.toc.currentItem().text() == "Step 4: Marker (ladder)"
    # tour: every step points at an existing part of the panel; the frame is removed again
    w.show_tour()
    t = w._tour
    css = {n: getattr(w, n).styleSheet() for _, _, n in help_dialog.TOUR if n}
    for k, (_, _, target) in enumerate(help_dialog.TOUR):
        t.go(k)
        if target:
            assert getattr(w, target) is t._framed and "E69F00" in getattr(w, target).styleSheet()
    t.close()
    assert all(getattr(w, n).styleSheet() == c for n, c in css.items())


def test_guide_is_plain_language():
    # user request: no technical jargon in the in-app help (the full manual keeps the details)
    from autogel.help_dialog import GUIDE_PATH, TOUR
    with open(GUIDE_PATH, encoding="utf-8") as f:
        text = f.read().lower() + " ".join(t + " " + x for t, x, _ in TOUR).lower()
    for word in ("deconvol", "rolling ball", "r²", "regression", "gaussian", "pipeline", "algorithm",
                 "signal/noise", "snr", "threshold", "parameter", "onnx", "pixel", "calibration curve",
                 "dynamic programming", "csv", "json"):
        assert word not in text, word

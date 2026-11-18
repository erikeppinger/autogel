"""napari editor for autogel.

Layers
------
gel        oriented image
saturated  red overlay of clipped pixels
lanes      rectangles, one per lane - move/resize/delete/add to edit lanes
bands      rectangles, one per band - the vertical extent is the integration range;
           draw a new rectangle inside a lane to add a band, delete to remove one

Every edit is applied automatically (short delay) and written to the log.
"""
from __future__ import annotations

import contextlib
import os
import warnings

import numpy as np
from qtpy.QtCore import QSettings, Qt, QTimer
from qtpy.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .analysis import GelAnalysis, Lane, Params
from .export import export_all, plot_lane
from .markers import DEFAULT_PRESET, PRESETS, parse_sizes, save_user_ladder

# Okabe & Ito (2008) colour-universal-design palette: distinguishable with
# protanopia, deuteranopia and tritanopia; "classic" kept for users who prefer it.
PALETTES = {
    "colour-blind safe (Okabe-Ito)": {"lane": "#56B4E9", "marker": "#E69F00", "band": "#F0E442",
                                      "user_band": "#CC79A7", "saturated": "#D55E00",
                                      "lane_text": "#FFFFFF", "band_text": "#F0E442"},
    "classic": {"lane": "cyan", "marker": "gold", "band": "lime", "user_band": "magenta",
                "saturated": "red", "lane_text": "yellow", "band_text": "lime"},
}
DEFAULT_PALETTE = "colour-blind safe (Okabe-Ito)"


# the profile dock can be collapsed to ~0 px height; matplotlib then warns on every redraw
warnings.filterwarnings("ignore", message="constrained_layout not applied")


def _rect(y0, x0, y1, x1):
    return np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]], float)


def _bbox(shape):
    s = np.asarray(shape)
    return (float(s[:, 0].min()), float(s[:, 1].min()), float(s[:, 0].max()), float(s[:, 1].max()))


class ProfileWidget(QWidget):
    """Matplotlib view of one lane profile with lane navigation and band hover.

    Callbacks (set by the owner):
      on_lane(i)       user stepped to lane i (buttons, spin box, wheel, arrow keys)
      on_hover(band)   mouse is over band index ``band`` (global index) or None
    """

    def __init__(self):
        super().__init__()
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        self.on_lane = lambda i: None
        self.on_hover = lambda b: None
        self.ga = None
        self.df = None
        self.lane = 0
        self._band_ids: list[int] = []
        self._hover: int | None = None
        self._hl = None

        bar = QHBoxLayout()
        self.btn_prev = QPushButton("◀")
        self.btn_next = QPushButton("▶")
        for btn in (self.btn_prev, self.btn_next):
            btn.setFixedWidth(32)
        self.sp_lane = QSpinBox()
        self.sp_lane.setRange(1, 1)
        self.sp_lane.setPrefix("lane ")
        self.lbl = QLabel("")
        hint = QLabel("wheel / ← → : change lane · hover a peak: highlight band")
        hint.setStyleSheet("color: gray")
        self.btn_prev.clicked.connect(lambda: self.step(-1))
        self.btn_next.clicked.connect(lambda: self.step(+1))
        self.sp_lane.valueChanged.connect(lambda v: self._goto(v - 1))
        bar.addWidget(self.btn_prev)
        bar.addWidget(self.sp_lane)
        bar.addWidget(self.btn_next)
        bar.addWidget(self.lbl)
        bar.addStretch()
        bar.addWidget(hint)

        self.fig = Figure(figsize=(6, 2.6), layout="constrained")  # no warnings when the dock is small
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.ax = self.fig.add_subplot(111)
        self.canvas.mpl_connect("motion_notify_event", self._on_move)
        self.canvas.mpl_connect("axes_leave_event", lambda e: self._set_hover(None))
        self.canvas.mpl_connect("figure_leave_event", lambda e: self._set_hover(None))
        self.canvas.mpl_connect("scroll_event", lambda e: self.step(-1 if e.button == "up" else +1))
        self.canvas.mpl_connect("key_press_event", self._on_key)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addLayout(bar)
        lay.addWidget(self.canvas)

    # ---------------------------------------------------------- navigation
    def step(self, d: int):
        if self.ga is not None and self.ga.lanes:
            self._goto(int(np.clip(self.lane + d, 0, len(self.ga.lanes) - 1)))

    def _goto(self, i: int):
        if self.ga is None or not self.ga.lanes or i == self.lane:
            return
        self.show_lane(self.ga, i, self.df)
        self.on_lane(i)

    def _on_key(self, e):
        if e.key in ("left", "up"):
            self.step(-1)
        elif e.key in ("right", "down"):
            self.step(+1)

    # --------------------------------------------------------------- hover
    def _on_move(self, e):
        if self.ga is None or e.inaxes is not self.ax or e.xdata is None:
            self._set_hover(None)
            return
        hit = None
        for k in self._band_ids:
            b = self.ga.bands[k]
            if b.y0 - 1 <= e.xdata <= b.y1 + 1:
                hit = k
                break
        self._set_hover(hit)

    def _set_hover(self, k):
        if k == self._hover:
            return
        self._hover = k
        if self._hl is not None:
            with contextlib.suppress(ValueError):   # artist already gone after ax.clear()
                self._hl.remove()
            self._hl = None
        if k is not None and self.ga is not None and k < len(self.ga.bands):
            b = self.ga.bands[k]
            self._hl = self.ax.axvspan(b.y0, b.y1, color="#CC79A7", alpha=0.35, zorder=0)
        self.canvas.draw_idle()
        self.on_hover(k)

    # ---------------------------------------------------------------- draw
    def show_lane(self, ga: GelAnalysis | None, i: int, df=None):
        self.ga, self.df = ga, df
        self._hover, self._hl = None, None
        self.ax.clear()
        n = len(ga.lanes) if ga is not None else 0
        self.lane = int(np.clip(i, 0, max(0, n - 1)))
        self.sp_lane.blockSignals(True)
        self.sp_lane.setRange(1, max(1, n))
        self.sp_lane.setValue(self.lane + 1)
        self.sp_lane.blockSignals(False)
        self.btn_prev.setEnabled(self.lane > 0)
        self.btn_next.setEnabled(self.lane < n - 1)
        self._band_ids = []
        if n:
            plot_lane(ga, self.lane, self.ax, df)
            self._band_ids = [k for k, b in enumerate(ga.bands) if b.lane == self.lane]
            self.lbl.setText(ga.lanes[self.lane].name)
        else:
            self.lbl.setText("")
        self.canvas.draw_idle()


class AutogelWidget(QWidget):
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        self.ga: GelAnalysis | None = None
        self.df = None
        self._updating = False
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(400)
        self._timer.timeout.connect(self._apply_edits)
        self._pending: set[str] = set()
        self._focus_lane: int | None = None
        self._build_ui()
        self.profile = ProfileWidget()
        self.profile.on_lane = self._profile_lane_changed
        self.profile.on_hover = self._profile_hover
        self.viewer.window.add_dock_widget(self.profile, name="lane profile", area="bottom")
        self._apply_mode()

    # ================================================================== UI
    def _build_ui(self):
        inner = QWidget()
        lay = QVBoxLayout(inner)
        self._expert_items: list = []   # (form_layout | None, widget) shown only in expert mode

        # --- mode switch (like "simple / expert" in slicer software)
        row = QHBoxLayout()
        row.addWidget(QLabel("Mode:"))
        self.btn_easy = QPushButton("Easy")
        self.btn_expert = QPushButton("Expert")
        self._mode_group = QButtonGroup(self)
        for btn in (self.btn_easy, self.btn_expert):
            btn.setCheckable(True)
            self._mode_group.addButton(btn)
            row.addWidget(btn)
        row.addStretch()
        lay.addLayout(row)
        mode = str(QSettings("autogel", "autogel").value("mode", "easy"))
        (self.btn_expert if mode == "expert" else self.btn_easy).setChecked(True)
        self.btn_easy.toggled.connect(lambda _: self._apply_mode(save=True))

        # --- file
        g = QGroupBox("1  Image")
        f = QVBoxLayout(g)
        row = QHBoxLayout()
        b = QPushButton("Open image…")
        b.clicked.connect(self._open_image_dialog)
        row.addWidget(b)
        b = QPushButton("Open session…")
        b.clicked.connect(self._open_session_dialog)
        row.addWidget(b)
        f.addLayout(row)
        self.lbl_file = QLabel("no image")
        self.lbl_file.setWordWrap(True)
        f.addWidget(self.lbl_file)
        self.lbl_warn = QLabel("")
        self.lbl_warn.setWordWrap(True)
        self.lbl_warn.setStyleSheet("color: orange")
        f.addWidget(self.lbl_warn)
        orient = QWidget()
        ol = QVBoxLayout(orient)
        ol.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        for txt, fn in (("Rotate 90°", self._rot90), ("Flip ↔", self._flip_lr),
                        ("Flip ↕", self._flip_ud)):
            b = QPushButton(txt)
            b.clicked.connect(fn)
            row.addWidget(b)
        ol.addLayout(row)
        row = QHBoxLayout()
        self.chk_invert = QCheckBox("dark bands (Coomassie)")
        self.chk_invert.toggled.connect(self._invert_toggled)
        row.addWidget(self.chk_invert)
        row.addWidget(QLabel("skew °"))
        self.sp_angle = QDoubleSpinBox()
        self.sp_angle.setRange(-10, 10)
        self.sp_angle.setSingleStep(0.25)
        self.sp_angle.editingFinished.connect(self._angle_changed)
        row.addWidget(self.sp_angle)
        ol.addLayout(row)
        f.addWidget(orient)
        self._expert_items.append((None, orient))
        lay.addWidget(g)

        # --- lanes
        g = QGroupBox("2  Lanes")
        f = QFormLayout(g)
        self.sp_lanes = QSpinBox()
        self.sp_lanes.setRange(0, 200)
        self.sp_lanes.setSpecialValueText("ask")
        self.sp_lanes.setToolTip("Number of lanes incl. marker. 'ask' = estimate and ask you to confirm.")
        f.addRow("number of lanes", self.sp_lanes)
        b = QPushButton("Detect lanes")
        b.clicked.connect(self.detect_lanes)
        f.addRow(b)
        self.ed_names = QLineEdit()
        self.ed_names.setPlaceholderText("lane names, comma separated (optional)")
        self.ed_names.editingFinished.connect(self._names_changed)
        f.addRow("names", self.ed_names)
        lay.addWidget(g)

        # --- bands
        g = QGroupBox("3  Bands")
        f = QFormLayout(g)
        self.cb_bg = QComboBox()
        self.cb_bg.addItems(["rolling_ball", "lane_edges", "none"])
        f.addRow("background", self.cb_bg)
        self.sp_radius = QDoubleSpinBox()
        self.sp_radius.setRange(0, 2000)
        self.sp_radius.setSpecialValueText("auto")
        f.addRow("rolling ball radius px", self.sp_radius)
        self.sp_snr = QDoubleSpinBox()
        self.sp_snr.setRange(1, 100)
        self.sp_snr.setValue(Params.min_snr)
        f.addRow("min. signal/noise", self.sp_snr)
        self.sp_cov = QDoubleSpinBox()
        self.sp_cov.setRange(0, 1)
        self.sp_cov.setSingleStep(0.05)
        self.sp_cov.setValue(Params.min_coverage)
        f.addRow("min. lane coverage", self.sp_cov)
        for w in (self.cb_bg, self.sp_radius, self.sp_snr, self.sp_cov):
            self._expert_items.append((f, w))
        b = QPushButton("Detect bands (all lanes)")
        b.clicked.connect(self.detect_bands)
        f.addRow(b)
        b = QPushButton("Re-detect bands in selected lane")
        b.clicked.connect(self._redetect_selected)
        f.addRow(b)
        self._expert_items.append((f, b))
        hint = QLabel("Edit in the viewer: select the 'lanes' or 'bands' layer, then move/resize/"
                      "delete rectangles, or draw new ones (rectangle tool). Band height = "
                      "integration range. Changes apply automatically and are logged.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: gray")
        f.addRow(hint)
        lay.addWidget(g)

        # --- marker
        g = QGroupBox("4  Marker / calibration")
        f = QFormLayout(g)
        self.cb_preset = QComboBox()
        self.cb_preset.addItems(list(PRESETS) + ["custom"])
        self.cb_preset.setCurrentText(DEFAULT_PRESET)
        self.cb_preset.currentTextChanged.connect(self._preset_changed)
        f.addRow("ladder", self.cb_preset)
        self.lbl_ladder = QLabel("")
        self.lbl_ladder.setWordWrap(True)
        self.lbl_ladder.setOpenExternalLinks(True)
        self.lbl_ladder.setStyleSheet("color: gray")
        f.addRow(self.lbl_ladder)
        self.ed_sizes = QLineEdit()
        f.addRow("sizes (top→bottom)", self.ed_sizes)
        self.ed_masses = QLineEdit()
        self.ed_masses.setPlaceholderText("optional: mass per band (datasheet)")
        f.addRow("masses", self.ed_masses)
        b = QPushButton("Save sizes/masses as my ladder…")
        b.clicked.connect(self._save_ladder)
        f.addRow(b)
        for w in (self.ed_sizes, self.ed_masses, b):
            self._expert_items.append((f, w))
        self.ed_marker_lanes = QLineEdit()
        self.ed_marker_lanes.setPlaceholderText("e.g. 1, 20")
        f.addRow("marker lane(s)", self.ed_marker_lanes)
        self.sp_load = QDoubleSpinBox()
        self.sp_load.setRange(0, 1e6)
        self.sp_load.setDecimals(3)
        self.sp_load.setSuffix(" µg")
        self.sp_load.setSpecialValueText("as datasheet")
        self.sp_load.setToolTip("Amount of ladder DNA you loaded (µg). µl → µg: volume × stock concentration")
        f.addRow("ladder loaded", self.sp_load)
        self.sp_ref = QDoubleSpinBox()
        self.sp_ref.setRange(0, 1e6)
        self.sp_ref.setDecimals(3)
        self.sp_ref.setSuffix(" µg")
        self.sp_ref.setSpecialValueText("from library")
        self.sp_ref.setToolTip("Ladder amount the datasheet masses refer to; filled from the ladder library")
        f.addRow("datasheet load", self.sp_ref)
        self.ed_loads = QLineEdit()
        self.ed_loads.setPlaceholderText("sample µl per lane, e.g. 2:5, 3:5 (→ ng/µl)")
        f.addRow("sample loads", self.ed_loads)
        for w in (self.sp_ref, self.ed_loads):
            self._expert_items.append((f, w))
        b = QPushButton("Calibrate && quantify")
        b.clicked.connect(self.quantify)
        f.addRow(b)
        self.lbl_cal = QLabel("")
        self.lbl_cal.setWordWrap(True)
        f.addRow(self.lbl_cal)
        lay.addWidget(g)
        self._preset_changed(DEFAULT_PRESET)

        # --- display
        g = QGroupBox("Display")
        f = QFormLayout(g)
        self.cb_palette = QComboBox()
        self.cb_palette.addItems(list(PALETTES))
        self.cb_palette.setCurrentText(DEFAULT_PALETTE)
        self.cb_palette.currentTextChanged.connect(self._palette_changed)
        f.addRow("colours", self.cb_palette)
        self.chk_all_labels = QCheckBox("label all bands (default: marker + selected lane)")
        self.chk_all_labels.toggled.connect(lambda _: self._update_band_labels())
        f.addRow(self.chk_all_labels)
        lay.addWidget(g)
        self._expert_items.append((None, g))

        # --- results
        g = QGroupBox("5  Results")
        f = QVBoxLayout(g)
        self.table = QTableWidget()
        self.table.setMinimumHeight(220)
        f.addWidget(self.table)
        row = QHBoxLayout()
        b = QPushButton("Export…")
        b.clicked.connect(self.export)
        row.addWidget(b)
        b = QPushButton("Save session…")
        b.clicked.connect(self.save_session)
        row.addWidget(b)
        f.addLayout(row)
        lay.addWidget(g)
        lay.addStretch()

        sc = QScrollArea()
        sc.setWidgetResizable(True)
        sc.setWidget(inner)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(sc)
        self.setMinimumWidth(380)

    def _apply_mode(self, save: bool = False):
        expert = self.btn_expert.isChecked()
        for form, w in self._expert_items:
            if form is None:
                w.setVisible(expert)
            else:
                form.setRowVisible(w, expert)
        if save:
            QSettings("autogel", "autogel").setValue("mode", "expert" if expert else "easy")

    # ============================================================ helpers
    def _params(self) -> Params:
        p = self.ga.params if self.ga else Params()
        p.bg_method = self.cb_bg.currentText()
        if self.sp_radius.value() > 0:
            p.bg_radius = self.sp_radius.value()
        p.min_snr = self.sp_snr.value()
        p.min_coverage = self.sp_cov.value()
        p.n_lanes = self.sp_lanes.value() or None
        p.invert = self.chk_invert.isChecked()
        return p

    def _info(self, msg):
        self.viewer.status = msg
        try:
            from napari.utils.notifications import show_info
            show_info(msg)
        except Exception:  # pragma: no cover
            pass

    def _need(self) -> bool:
        if self.ga is None:
            QMessageBox.information(self, "autogel", "Open an image first.")
            return False
        return True

    # ============================================================ loading
    def _open_image_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open gel image", "",
                                              "Images (*.tif *.tiff *.png *.jpg *.jpeg *.bmp);;All (*)")
        if path:
            self.open_image(path)

    def _open_session_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open session", "", "Session (*.json)")
        if not path:
            return
        try:
            ga = GelAnalysis.load_session(path)
        except FileNotFoundError:
            img, _ = QFileDialog.getOpenFileName(self, "Image of this session not found - locate it")
            if not img:
                return
            ga = GelAnalysis.load_session(path, img)
        self._set_analysis(ga)
        self._sync_controls_from_params()

    def open_image(self, path: str, orientation: dict | None = None):
        p = self._params() if self.ga else Params(invert=None)
        p.n_lanes = self.sp_lanes.value() or None
        if not self.ga:
            p.invert = None
        self.ga = None
        ga = GelAnalysis(path, p, orientation=orientation)
        self._set_analysis(ga)

    def _set_analysis(self, ga: GelAnalysis):
        self.ga = ga
        self.df = None
        self._focus_lane = None
        self._updating = True
        try:
            self.chk_invert.blockSignals(True)
            self.chk_invert.setChecked(bool(ga.gel.orientation.get("invert")))
            self.chk_invert.blockSignals(False)
            self.sp_angle.setValue(ga.gel.orientation.get("angle", 0.0))
            for name in ("gel", "saturated", "highlight", "lanes", "bands"):
                if name in self.viewer.layers:
                    self.viewer.layers.remove(name)
            a = ga.data
            lo, hi = np.percentile(a, [0.5, 99.7])
            self.viewer.add_image(a, name="gel", colormap="gray", contrast_limits=(lo, hi))
            self.viewer.add_image(ga.gel.raw_max.astype(np.float32), name="saturated",
                                  colormap=self._sat_colormap(), blending="translucent",
                                  contrast_limits=(0, 1), opacity=1.0)
            self.viewer.add_shapes(name="highlight", shape_type="rectangle", edge_color="#CC79A7",
                                   face_color="transparent", edge_width=6, opacity=1.0)
            lanes = self.viewer.add_shapes(name="lanes", shape_type="rectangle",
                                           edge_color=self.colors["lane"],
                                           face_color="transparent", edge_width=2)
            bands = self.viewer.add_shapes(name="bands", shape_type="rectangle",
                                           edge_color=self.colors["band"],
                                           face_color="transparent", edge_width=1.5)
            lanes.events.data.connect(lambda e: self._schedule("lanes"))
            bands.events.data.connect(lambda e: self._schedule("bands"))
            lanes.events.highlight.connect(lambda e: self._lane_selected())
            bands.events.highlight.connect(lambda e: self._band_selected())
            self.viewer.reset_view()
        finally:
            self._updating = False
        self.lbl_file.setText(os.path.basename(ga.path) + f"  ({a.shape[1]}×{a.shape[0]} px, {ga.gel.dtype})")
        self.lbl_warn.setText("\n".join("⚠ " + w for w in ga.gel.warnings))
        self._refresh_layers()

    def _sync_controls_from_params(self):
        p = self.ga.params
        self.sp_lanes.setValue(p.n_lanes or 0)
        self.cb_bg.setCurrentText(p.bg_method)
        self.sp_snr.setValue(p.min_snr)
        self.sp_cov.setValue(p.min_coverage)
        if p.marker_preset in PRESETS:
            self.cb_preset.setCurrentText(p.marker_preset)
        if p.marker_sizes:
            self.ed_sizes.setText(", ".join(f"{s:g}" for s in p.marker_sizes))
        self.ed_masses.setText(", ".join(f"{m:g}" for m in p.marker_masses) if p.marker_masses else "")
        self.ed_marker_lanes.setText(", ".join(map(str, p.marker_lanes)))
        self.ed_names.setText(", ".join(l.name for l in self.ga.lanes))

    # ------------------------------------------------------- orientation
    def _reorient(self, change):
        if not self._need():
            return
        o = dict(self.ga.gel.orientation)
        change(o)
        self.ga._log("orientation_change", o)
        log = self.ga.log
        self.open_image(self.ga.path, orientation=o)
        self.ga.log = log + self.ga.log[1:]

    def _rot90(self):
        self._reorient(lambda o: o.update(rot90_k=(o.get("rot90_k", 0) + 1) % 4))

    def _flip_lr(self):
        self._reorient(lambda o: o.update(flip_lr=not o.get("flip_lr", False)))

    def _flip_ud(self):
        self._reorient(lambda o: o.update(flip_ud=not o.get("flip_ud", False)))

    def _angle_changed(self):
        if self.ga and abs(self.sp_angle.value() - self.ga.gel.orientation.get("angle", 0)) > 1e-6:
            self._reorient(lambda o: o.update(angle=self.sp_angle.value()))

    def _invert_toggled(self, on):
        if self.ga and bool(self.ga.gel.orientation.get("invert")) != on:
            self._reorient(lambda o: o.update(invert=on))

    # ------------------------------------------------------------- lanes
    def detect_lanes(self):
        if not self._need():
            return
        p = self._params()
        asked = False
        if not p.n_lanes:
            # estimate first, then let the user confirm the real number (internal check)
            p.n_lanes = None
            est = len(self.ga.detect_lanes())
            n, ok = QInputDialog.getInt(self, "Number of lanes",
                                        f"Detected {est} lanes.\nHow many lanes (incl. marker) are on the gel?",
                                        est, 1, 200)
            if not ok:
                return
            self.sp_lanes.setValue(n)
            p.n_lanes = n
            asked = True
            if n == est:
                self.ga._log("confirm_lanes", {"n_lanes": n}, "user confirmed automatic lane count")
                self._refresh_layers()
                self._names_changed()
                return
        self.ga.detect_lanes()
        found = len(self.ga.lanes)
        if found != p.n_lanes:
            QMessageBox.warning(self, "autogel",
                                f"Expected {p.n_lanes} lanes but could only place {found}. "
                                "Please add/adjust lanes in the 'lanes' layer.")
        if asked:
            self.ga._log("confirm_lanes", {"n_lanes": p.n_lanes}, "user corrected lane count")
        self._refresh_layers()
        self._names_changed()

    def _names_changed(self):
        if not self.ga or not self.ed_names.text().strip():
            return
        names = [n.strip() for n in self.ed_names.text().split(",")]
        names += [str(i + 1) for i in range(len(names), len(self.ga.lanes))]
        self.ga.set_lane_names(names[: len(self.ga.lanes)])
        self._refresh_layers()

    # ------------------------------------------------------------- bands
    def detect_bands(self):
        if not self._need():
            return
        if not self.ga.lanes:
            self.detect_lanes()
        self._params()
        self.ga._profiles = {}
        self.ga.detect_bands()
        self._refresh_layers()
        self.quantify(silent=True)

    def _selected_lane(self) -> int | None:
        if "lanes" in self.viewer.layers:
            sel = list(self.viewer.layers["lanes"].selected_data)
            if sel:
                return sel[0]
        return None

    def _redetect_selected(self):
        if not self._need():
            return
        i = self._selected_lane()
        if i is None:
            QMessageBox.information(self, "autogel", "Select a lane in the 'lanes' layer first.")
            return
        self._params()
        self.ga._profiles = {}
        self.ga.detect_bands([i])
        self._refresh_layers()
        self.quantify(silent=True)

    # ------------------------------------------------------------- edits
    def _schedule(self, which):
        if self._updating or self.ga is None:
            return
        self._pending.add(which)
        self._timer.start()

    def _apply_edits(self):
        if self.ga is None:
            return
        pending, self._pending = self._pending, set()
        if "lanes" in pending:
            shapes = self.viewer.layers["lanes"].data
            old = {(round(l.x), round(l.width), round(l.y0), round(l.y1)): l for l in self.ga.lanes}
            lanes = []
            for s in shapes:
                y0, x0, y1, x1 = _bbox(s)
                key = (round((x0 + x1) / 2), round(x1 - x0), round(y0), round(y1))
                prev = old.get(key)
                lanes.append(Lane(x=(x0 + x1) / 2, width=x1 - x0, y0=y0, y1=y1,
                                  name=prev.name if prev else "",
                                  source=prev.source if prev else "user"))
            self.ga.set_lanes(lanes, reason="edited in viewer")
        if "bands" in pending and "lanes" not in pending:
            old = {(b.lane, round(b.y0), round(b.y1)): b for b in self.ga.bands}
            bands, orphan = [], 0
            for s in self.viewer.layers["bands"].data:
                y0, x0, y1, x1 = _bbox(s)
                y1 -= 1  # rectangles are drawn to b.y1 + 1 (inclusive last row)
                xc = (x0 + x1) / 2
                hit = [i for i, l in enumerate(self.ga.lanes) if l.x0 <= xc <= l.x1]
                if not hit:
                    orphan += 1
                    continue
                i = hit[0]
                prev = old.get((i, round(y0), round(y1)))
                bands.append(prev if prev else self.ga.band_for(i, y0, y1, source="user"))
            self.ga.set_bands(bands, reason="edited in viewer")
            if orphan:
                self._info(f"{orphan} band rectangle(s) outside any lane were ignored")
        self._refresh_layers()
        self.quantify(silent=True)

    # ------------------------------------------------------------ display
    @property
    def colors(self) -> dict:
        return PALETTES[self.cb_palette.currentText()]

    def _sat_colormap(self):
        from napari.utils.colormaps import Colormap
        return Colormap([[0, 0, 0, 0], self.colors["saturated"]], name="autogel_sat")

    def _palette_changed(self, _=None):
        if self.ga is None:
            return
        if "saturated" in self.viewer.layers:
            self.viewer.layers["saturated"].colormap = self._sat_colormap()
        self._refresh_layers()

    def _refresh_layers(self):
        if self.ga is None:
            return
        self._updating = True
        col = self.colors
        try:
            ga = self.ga
            marker = {L - 1 for L in ga.params.marker_lanes}
            lanes = self.viewer.layers["lanes"]
            lanes.data = []
            if ga.lanes:
                lanes.add_rectangles([_rect(l.y0, l.x0, l.y1, l.x1) for l in ga.lanes],
                                     edge_color=[col["marker"] if i in marker else col["lane"]
                                                 for i in range(len(ga.lanes))])
                lanes.features = {"name": [l.name for l in ga.lanes]}
                # lane names centred ABOVE the lanes so they never sit on band labels
                h = float(np.median([l.y1 - l.y0 for l in ga.lanes]))
                lanes.text = {"string": "{name}", "anchor": "center", "color": col["lane_text"],
                              "size": 8, "translation": [-(h / 2 + 25), 0]}
            bands = self.viewer.layers["bands"]
            bands.data = []
            if ga.bands:
                bands.add_rectangles([_rect(b.y0, ga.lanes[b.lane].x0, b.y1 + 1, ga.lanes[b.lane].x1)
                                      for b in ga.bands],
                                     edge_color=[col["band"] if b.source == "auto" else col["user_band"]
                                                 for b in ga.bands])
        finally:
            self._updating = False
        self._update_band_labels()
        self._lane_selected()

    def _update_band_labels(self):
        """Size labels only where readable: marker lanes and the selected lane
        (or every band when 'label all bands' is ticked)."""
        if self.ga is None or "bands" not in self.viewer.layers:
            return
        ga = self.ga
        bands = self.viewer.layers["bands"]
        if not ga.bands or len(bands.data) != len(ga.bands):
            return
        sizes = {}
        if self.df is not None and len(self.df):
            sc = [c for c in self.df.columns if c.startswith("size_")]
            if sc:
                sizes = {(r["lane"] - 1, r["y_px"]): r[sc[0]] for _, r in self.df.iterrows()}
        show = {L - 1 for L in ga.params.marker_lanes}
        if self._focus_lane is not None:
            show.add(self._focus_lane)
        labels = []
        for b in ga.bands:
            v = sizes.get((b.lane, round(b.y, 1)))
            ok = self.chk_all_labels.isChecked() or b.lane in show
            labels.append(f"{v:.0f}" if ok and v is not None and np.isfinite(v) else "")
        self._updating = True
        try:
            bands.features = {"label": labels}
            bands.text = {"string": "{label}", "anchor": "center", "color": self.colors["band_text"],
                          "size": 7}
        finally:
            self._updating = False

    def _lane_selected(self):
        i = self._selected_lane()
        if i is not None and i != self._focus_lane:
            self._focus_lane = i
            self._update_band_labels()
        lane = self._focus_lane if self._focus_lane is not None else 0
        self.profile.show_lane(self.ga, lane, self.df)
        self._highlight(lane, None)

    def _band_selected(self):
        if self.ga is None or "bands" not in self.viewer.layers:
            return
        sel = list(self.viewer.layers["bands"].selected_data)
        if sel and sel[0] < len(self.ga.bands):
            lane = self.ga.bands[sel[0]].lane
            if lane != self._focus_lane:
                self._focus_lane = lane
                self._update_band_labels()
                self.profile.show_lane(self.ga, lane, self.df)
            self._highlight(lane, sel[0])

    def _profile_lane_changed(self, i: int):
        """Lane stepped in the profile panel (buttons, wheel, keys)."""
        self._focus_lane = i
        self._update_band_labels()
        self._highlight(i, None)

    def _profile_hover(self, band: int | None):
        lane = self.profile.lane if self.ga and self.ga.lanes else None
        self._highlight(lane, band)

    def _highlight(self, lane: int | None, band: int | None):
        """Outline the current lane (and hovered band) in the image; never touches edits."""
        if self.ga is None or "highlight" not in self.viewer.layers:
            return
        hl = self.viewer.layers["highlight"]
        shapes, widths = [], []
        if lane is not None and 0 <= lane < len(self.ga.lanes):
            l = self.ga.lanes[lane]
            shapes.append(_rect(l.y0, l.x0 - 3, l.y1, l.x1 + 3))
            widths.append(3)
        if band is not None and 0 <= band < len(self.ga.bands):
            b = self.ga.bands[band]
            l = self.ga.lanes[b.lane]
            shapes.append(_rect(b.y0 - 3, l.x0 - 6, b.y1 + 4, l.x1 + 6))
            widths.append(7)
        self._updating = True
        try:
            hl.data = []
            if shapes:
                hl.add_rectangles(shapes, edge_width=widths, edge_color="#CC79A7",
                                  face_color="transparent")
        finally:
            self._updating = False

    # ------------------------------------------------------- calibration
    def _preset_changed(self, name):
        if name not in PRESETS:
            self.lbl_ladder.setText("custom ladder: enter sizes (and masses) below")
            return
        pre = PRESETS[name]
        self.ed_sizes.setText(", ".join(f"{s:g}" for s in pre["sizes"]))
        self.ed_masses.setText(", ".join(f"{m:g}" for m in pre["masses"]) if pre.get("masses") else "")
        ref = pre.get("reference_load") or {}
        self.sp_ref.setValue(float(ref.get("amount") or 0))
        bits = []
        if pre.get("user_defined"):
            bits.append("your ladder")
        elif pre.get("verified"):
            bits.append(f"✔ checked against vendor data {pre.get('checked', '')}")
        elif pre.get("source_db") == "ApE":
            bits.append("⚠ imported from ApE (sizes only) - compare with your datasheet")
        else:
            bits.append("⚠ not yet verified - compare with your datasheet")
        if ref.get("text"):
            bits.append(f"masses: {ref['text']}")
        elif not pre.get("masses"):
            bits.append("no band masses → no ng estimate")
        if pre.get("note"):
            bits.append(pre["note"])
        if pre.get("source"):
            bits.append(f'<a href="{pre["source"]}">source</a>')
        self.lbl_ladder.setText("<br>".join(bits))

    def _save_ladder(self):
        name, ok = QInputDialog.getText(self, "Save ladder", "Name for this ladder:")
        if not ok or not name.strip():
            return
        unit, ok = QInputDialog.getItem(self, "Save ladder", "Unit:", ["bp", "kDa"], 0, False)
        if not ok:
            return
        try:
            sizes = parse_sizes(self.ed_sizes.text())
            masses = [float(t) for t in self.ed_masses.text().replace(";", ",").split(",") if t.strip()] or None
            path = save_user_ladder(name.strip(), sizes, unit, masses)
        except ValueError as e:
            QMessageBox.warning(self, "autogel", str(e))
            return
        self.cb_preset.blockSignals(True)
        self.cb_preset.clear()
        self.cb_preset.addItems(list(PRESETS) + ["custom"])
        self.cb_preset.setCurrentText(name.strip())
        self.cb_preset.blockSignals(False)
        self._info(f"ladder saved to {path}")

    def _apply_marker(self):
        txt = self.ed_marker_lanes.text().strip()
        if not txt:
            return False
        try:
            lanes = [int(t) for t in txt.replace(";", ",").replace(" ", ",").split(",") if t.strip()]
            sizes = parse_sizes(self.ed_sizes.text())
            masses = [float(t) for t in self.ed_masses.text().replace(";", ",").split(",") if t.strip()] or None
        except ValueError as e:
            raise ValueError(f"Could not read marker lanes / sizes / masses: {e}") from e
        bad = [L for L in lanes if not 1 <= L <= len(self.ga.lanes)]
        if bad:
            self._info(f"marker lane(s) {bad} do not exist (gel has {len(self.ga.lanes)} lanes)")
        if masses is not None and len(masses) != len(sizes):
            # non-modal: this runs after every edit
            self._info(f"{len(sizes)} sizes but {len(masses)} masses - masses ignored")
            masses = None
        preset = self.cb_preset.currentText()
        unit = PRESETS[preset]["unit"] if preset in PRESETS else self.ga.params.marker_unit
        self.ga.set_marker(lanes, preset=preset, sizes=sizes, masses=masses, unit=unit,
                           load=self.sp_load.value() or None, ref_load=self.sp_ref.value() or None)
        loads, bad = {}, []
        for part in self.ed_loads.text().split(","):
            if not part.strip():
                continue
            k, _, v = part.partition(":")
            try:
                loads[str(int(k.strip()))] = float(v)
            except ValueError:
                bad.append(part.strip())
        if bad:
            self._info(f"ignored sample loads (use lane:µl): {', '.join(bad)}")
        self.ga.params.sample_loads_ul = loads
        return True

    def quantify(self, silent=False):
        if not self._need():
            return
        try:
            if not silent or self.ga.params.marker_lanes:
                self._apply_marker()
            self.df = self.ga.quantify()
        except Exception as e:  # show, don't crash the viewer
            QMessageBox.critical(self, "autogel", f"Quantification failed:\n{e}")
            return
        msg = []
        for c in self.ga.calibrations:
            msg += c.warnings
        if self.ga.mass_fit:
            mf = self.ga.mass_fit
            msg.append(f"mass curve: R² = {mf['r2']:.3f} (n={mf['n']}, through origin)")
        self.lbl_cal.setText("\n".join(msg))
        self._fill_table()
        self._refresh_layers()

    def _fill_table(self):
        df = self.df
        self.table.clear()
        if df is None or not len(df):
            self.table.setRowCount(0)
            return
        cols = [c for c in df.columns if c not in ("y0_px", "y1_px", "raw_volume", "background")]
        self.table.setColumnCount(len(cols))
        self.table.setRowCount(len(df))
        self.table.setHorizontalHeaderLabels(cols)
        for r, (_, row) in enumerate(df[cols].iterrows()):
            for c, v in enumerate(row):
                it = QTableWidgetItem("" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
                it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if row.get("saturated"):
                    it.setForeground(Qt.GlobalColor.red)
                self.table.setItem(r, c, it)
        self.table.resizeColumnsToContents()

    # ------------------------------------------------------------ output
    def export(self):
        if not self._need():
            return
        d = QFileDialog.getExistingDirectory(self, "Export folder", os.path.dirname(self.ga.path))
        if not d:
            return
        try:
            self._apply_marker()
            out = export_all(self.ga, d)
        except (ValueError, OSError) as e:   # bad input or file locked (e.g. xlsx open in Excel)
            QMessageBox.critical(self, "autogel", f"Export failed:\n{e}")
            return
        self._info(f"exported to {out}")

    def save_session(self):
        if not self._need():
            return
        stem = os.path.splitext(os.path.basename(self.ga.path))[0]
        path, _ = QFileDialog.getSaveFileName(self, "Save session", stem + "_session.json", "Session (*.json)")
        if path:
            self.ga.save_session(path)
            self._info(f"session saved: {path}")


def launch(path: str | None = None):
    import napari
    from napari.qt import get_qapp
    get_qapp()  # napari >= 0.10 creates the QApplication lazily; widgets need it first
    viewer = napari.Viewer(title="autogel")
    w = AutogelWidget(viewer)
    viewer.window.add_dock_widget(w, name="autogel", area="right")
    if path:
        w.open_image(path)
    napari.run()
    return viewer

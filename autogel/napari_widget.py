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
    QMenu,
    QMessageBox,
    QPlainTextEdit,
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
from .help_dialog import MANUAL_URL, REPORT_URL, HelpDialog, TourDialog, open_url
from .labels import declutter, plan_labels
from .markers import DEFAULT_PRESET, PRESETS, RESOLUTION_GUIDE, parse_sizes, suggest_ladders
from .report import describe

# Okabe & Ito (2008) colour-universal-design palette: distinguishable with
# protanopia, deuteranopia and tritanopia; "classic" kept for users who prefer it.
PALETTES = {
    "colour-blind safe (Okabe-Ito)": {"lane": "#56B4E9", "marker": "#E69F00", "band": "#F0E442",
                                      "user_band": "#CC79A7", "saturated": "#D55E00",
                                      "lane_text": "#FFFFFF", "band_text": "#F0E442", "label_text": "#FFFFFF"},
    "classic": {"lane": "cyan", "marker": "gold", "band": "lime", "user_band": "magenta",
                "saturated": "red", "lane_text": "yellow", "band_text": "lime", "label_text": "white"},
}
DEFAULT_PALETTE = "colour-blind safe (Okabe-Ito)"
MASS_MODELS = {"size-dependent": "size", "one slope through origin": "linear"}
AI_BAND = "#009E73"        # bands added by the AI second opinion (Okabe-Ito bluish green)
DOUBT_BAND = "#999999"     # doubtful bands: unlikely (band checker) or not confirmed by the AI
DOUBT_PROBABILITY = 0.1
NORM_METHODS = {"none": "none", "loading-control band": "loading_control", "total lane signal": "total_lane"}


# the profile dock can be collapsed to ~0 px height; matplotlib then warns on every redraw
warnings.filterwarnings("ignore", message="constrained_layout not applied")


def _rect(y0, x0, y1, x1):
    return np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]], float)


def _lane_shape(l: Lane) -> np.ndarray:
    """What is drawn for a lane: a rectangle, or a polygon following a curved lane."""
    return l.outline() if l.bend else _rect(l.y0, l.x0, l.y1, l.x1)


def _lane_from_polygon(v: np.ndarray) -> Lane | None:
    """Inverse of Lane.outline(): left edge top->bottom, right edge bottom->top."""
    v = np.asarray(v, float)
    if len(v) < 4 or len(v) % 2:
        return None
    n = len(v) // 2
    left, right = v[:n], v[n:][::-1]
    ys = (left[:, 0] + right[:, 0]) / 2
    cs = (left[:, 1] + right[:, 1]) / 2
    width = float(np.median(right[:, 1] - left[:, 1]))
    if width <= 0 or np.any(np.diff(ys) <= 0):
        return None
    x = float(cs[0])
    return Lane(x=x, width=width, y0=float(ys[0]), y1=float(ys[-1]),
                bend=[[float(y), float(c - x)] for y, c in zip(ys, cs)])


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
        # Help sits next to the mode buttons: the right edge of a narrow dock can be cut off
        self.btn_help = QPushButton("Help ▾")
        self.btn_help.setToolTip("Quick tour, user guide (F1), full manual, report a problem")
        m = QMenu(self.btn_help)
        m.addAction("Quick tour", self.show_tour)
        m.addAction("User guide  (F1)", self.show_guide)
        m.addSeparator()
        m.addAction("Full manual (online)", lambda: open_url(MANUAL_URL))
        m.addAction("Report a problem or idea (online)", lambda: open_url(REPORT_URL))
        self.btn_help.setMenu(m)
        row.addWidget(self.btn_help)
        row.addStretch()
        lay.addLayout(row)
        from qtpy.QtGui import QKeySequence, QShortcut
        QShortcut(QKeySequence("F1"), self, activated=self.show_guide)
        # live activity line (always visible) + full step log (toggle)
        row = QHBoxLayout()
        self.lbl_status = QLabel("Open a gel image to start. New here? Help ▸ Quick tour.")
        self.lbl_status.setWordWrap(True)
        self.lbl_status.setStyleSheet("color: #56B4E9")
        row.addWidget(self.lbl_status, 1)
        self.btn_log = QPushButton("Log ▸")
        self.btn_log.setCheckable(True)
        self.btn_log.setFixedWidth(64)
        row.addWidget(self.btn_log)
        lay.addLayout(row)
        self.txt_log = QPlainTextEdit()
        self.txt_log.setReadOnly(True)
        self.txt_log.setMaximumBlockCount(2000)
        self.txt_log.setMinimumHeight(140)
        self.txt_log.setVisible(False)
        lay.addWidget(self.txt_log)
        self.btn_log.toggled.connect(self._toggle_log)
        mode = str(QSettings("autogel", "autogel").value("mode", "easy"))
        (self.btn_expert if mode == "expert" else self.btn_easy).setChecked(True)
        self.btn_easy.toggled.connect(lambda _: self._apply_mode(save=True))

        # --- file
        g = QGroupBox("1  Image")
        self.grp_image = g
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
        row = QHBoxLayout()
        row.addWidget(QLabel("channel"))
        self.cb_channel = QComboBox()
        self.cb_channel.setToolTip("Which image is analysed: colour channels, pages of a multi-page TIFF, "
                                   "added images, or 'all channels' to detect the bands of every colour.")
        self.cb_channel.currentIndexChanged.connect(self._channel_changed)
        row.addWidget(self.cb_channel, 1)
        b = QPushButton("Add image…")
        b.setToolTip("Add another image of the same gel/membrane (same camera position), e.g. the "
                     "white-light photo showing a prestained marker, or a total-protein stain.")
        b.clicked.connect(self._add_image_dialog)
        row.addWidget(b)
        f.addLayout(row)
        # rotate / flip stay visible in Easy mode: a wrong automatic orientation must be fixable
        row = QHBoxLayout()
        for txt, fn, tip in (("Rotate 90°", self._rot90, "Lanes must run top to bottom"),
                             ("Flip ↔", self._flip_lr, "Mirror left/right"),
                             ("Flip ↕", self._flip_ud, "Wells must be at the top")):
            b = QPushButton(txt)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            row.addWidget(b)
        f.addLayout(row)
        orient = QWidget()
        ol = QVBoxLayout(orient)
        ol.setContentsMargins(0, 0, 0, 0)
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
        self.grp_lanes = g
        f = QFormLayout(g)
        self.sp_lanes = QSpinBox()
        self.sp_lanes.setRange(0, 200)
        self.sp_lanes.setSpecialValueText("ask")
        self.sp_lanes.setToolTip("Number of lanes incl. marker. 'ask' = estimate and ask you to confirm.")
        f.addRow("number of lanes", self.sp_lanes)
        self.sp_tiers = QSpinBox()
        self.sp_tiers.setRange(1, 4)
        self.sp_tiers.setToolTip("Gels stacked in one image (one comb per gel row). The split line "
                                 "is found automatically and can be dragged in the 'gel rows' layer.")
        self.sp_tiers.valueChanged.connect(self._tiers_changed)
        f.addRow("gel rows (combs)", self.sp_tiers)
        self.sp_cols = QSpinBox()
        self.sp_cols.setRange(1, 4)
        self.sp_cols.setToolTip("Gels lying side by side in one image. The split is found automatically "
                                "and can be dragged in the 'gel splits' layer.")
        self.sp_cols.valueChanged.connect(self._cols_changed)
        f.addRow("gels side by side", self.sp_cols)
        self.chk_curved = QCheckBox("follow curved lanes")
        self.chk_curved.setChecked(Params.curved_lanes)
        self.chk_curved.setToolTip("Trace each lane down the gel and follow it when it clearly bends "
                                   "(fanned-out lanes). Untick for strictly straight lanes.")
        f.addRow(self.chk_curved)
        self._expert_items.append((f, self.chk_curved))
        self.chk_crop = QCheckBox("only the gel (auto-crop)")
        self.chk_crop.setChecked(Params.auto_crop)
        self.chk_crop.setToolTip("Ignore the tray, film edge or frame around the gel / membrane. The analysed "
                                 "region is shown as a rectangle (layer 'analysed region'); drag or resize it "
                                 "to set it by hand.")
        f.addRow(self.chk_crop)
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
        self.grp_bands = g
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
        self.sp_gelrel = QDoubleSpinBox()
        self.sp_gelrel.setRange(0, 0.5)
        self.sp_gelrel.setSingleStep(0.01)
        self.sp_gelrel.setValue(Params.min_gel_rel)
        self.sp_gelrel.setSpecialValueText("off")
        self.sp_gelrel.setToolTip("Ignore bands weaker than this fraction of the gel's strong bands. "
                                  "About 0.03 removes noise 'bands' in empty lanes of noisy or JPEG images. "
                                  "Leave off if faint ladder bands sit next to very bright samples.")
        f.addRow("min. height vs. gel", self.sp_gelrel)
        self.chk_floor = QCheckBox("automatic noise floor")
        self.chk_floor.setChecked(Params.noise_floor)
        self.chk_floor.setToolTip("Drop bands that are both weak (< 7 % of the gel's strong bands) and close to "
                                  "the gel's noise (< 10x). Keeps faint ladder bands, removes noise 'bands' in "
                                  "empty lanes of JPEG/noisy images.")
        f.addRow(self.chk_floor)
        self.chk_deconv = QCheckBox("split touching bands by peak fit")
        self.chk_deconv.setChecked(Params.deconvolve)
        self.chk_deconv.setToolTip("Bands whose ranges touch share their signal by fitted Gaussian peak areas "
                                   "instead of being cut at the valley (more accurate for a faint band next to "
                                   "a strong one). Column 'deconvolved' in the table.")
        f.addRow(self.chk_deconv)
        self.chk_shoulder = QCheckBox("split shoulders (unresolved doublets)")
        self.chk_shoulder.setChecked(Params.split_shoulders)
        self.chk_shoulder.setToolTip("Split a band in two when it clearly has a shoulder. Off by default: "
                                     "slightly lopsided real bands are occasionally split too - check the result.")
        f.addRow(self.chk_shoulder)
        self.sp_checker = QDoubleSpinBox()
        self.sp_checker.setRange(0, 0.9)
        self.sp_checker.setSingleStep(0.05)
        self.sp_checker.setValue(Params.checker_filter)
        self.sp_checker.setSpecialValueText("off")
        self.sp_checker.setToolTip("Band checker (learned from hand-labelled gels): remove bands whose probability "
                                   "of being real is below this value. Never in marker lanes. 0.1 removed a third "
                                   "of the false bands on unseen gels while keeping 98 % of the real ones.")
        f.addRow("remove unlikely bands", self.sp_checker)
        self._expert_items.append((f, self.sp_checker))
        for w in (self.cb_bg, self.sp_radius, self.sp_snr, self.sp_cov, self.sp_gelrel, self.chk_floor,
                  self.chk_deconv, self.chk_shoulder):
            self._expert_items.append((f, w))
        b = QPushButton("Detect bands (all lanes)")
        b.clicked.connect(self.detect_bands)
        f.addRow(b)
        b = QPushButton("AI second opinion")
        b.setToolTip("GelGenie's band-recognition network (optional extra): bands it does not confirm turn grey, "
                     "bands it sees where autogel has none are added (green). Check and correct them as usual.")
        b.clicked.connect(self.ai_second_opinion)
        f.addRow(b)
        b = QPushButton("Next band to check  ▸")
        b.setToolTip("Jumps to the next band worth a look: grey bands autogel is unsure about (most "
                     "doubtful first), then bands the AI suggested (green). Keep it, move it, or press "
                     "Delete to remove it, then click again.")
        b.clicked.connect(lambda: self.next_band_to_check())
        self.btn_next_check = b
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
        self.grp_marker = g
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
        for w in (self.ed_sizes, self.ed_masses):
            self._expert_items.append((f, w))
        lad = QWidget()
        ll = QHBoxLayout(lad)
        ll.setContentsMargins(0, 0, 0, 0)
        b = QPushButton("Ladders…")
        b.setToolTip("Add, edit, import/export or share ladders; update the library")
        b.clicked.connect(self._open_ladders)
        ll.addWidget(b)
        b = QPushButton("Suggest ladder")
        b.setToolTip("Which library ladder matches the bands of the marker lane?")
        b.clicked.connect(self._suggest_ladder)
        ll.addWidget(b)
        f.addRow(lad)
        # valid size range (resolution range of the gel): automatic, by gel %, or custom
        self.cb_range = QComboBox()
        self.cb_range.addItems(["automatic (whole ladder)", *RESOLUTION_GUIDE, "custom"])
        self.cb_range.setToolTip("Only ladder bands inside this range are used for the size and mass "
                                 "curves; sample bands outside are flagged. Gel % values are typical "
                                 "guidelines - adjust to your system.")
        f.addRow("valid size range", self.cb_range)
        rng = QWidget()
        rl = QHBoxLayout(rng)
        rl.setContentsMargins(0, 0, 0, 0)
        self.sp_rmin = QDoubleSpinBox()
        self.sp_rmax = QDoubleSpinBox()
        for sp in (self.sp_rmin, self.sp_rmax):
            sp.setRange(0, 1e7)
            sp.setDecimals(0)
            sp.setSpecialValueText("–")
            sp.setKeyboardTracking(False)
        rl.addWidget(self.sp_rmin)
        rl.addWidget(QLabel("to"))
        rl.addWidget(self.sp_rmax)
        f.addRow("  min / max", rng)
        self.cb_range.currentTextChanged.connect(self._range_preset_changed)
        self.sp_rmin.valueChanged.connect(self._range_values_changed)
        self.sp_rmax.valueChanged.connect(self._range_values_changed)
        self.ed_marker_lanes = QLineEdit()
        self.ed_marker_lanes.setPlaceholderText("e.g. 1, 20")
        f.addRow("marker lane(s)", self.ed_marker_lanes)
        self.cb_marker_ch = QComboBox()
        self.cb_marker_ch.setToolTip("Read the marker lanes from another image, e.g. the white-light photo of a "
                                     "Western blot (the chemiluminescence image usually doesn't show the marker).")
        self.cb_marker_ch.currentIndexChanged.connect(self._marker_channel_changed)
        f.addRow("marker lanes from", self.cb_marker_ch)
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
        self.cb_mass = QComboBox()
        self.cb_mass.addItems(list(MASS_MODELS))
        self.cb_mass.setToolTip("How band intensity is turned into ng. 'size-dependent' learns from the ladder "
                                "that small (and very large) fragments give less signal per ng.")
        self.cb_mass.currentTextChanged.connect(lambda _t: self.ga and self.ga.params.marker_lanes
                                                and self.quantify(silent=True))
        f.addRow("mass curve", self.cb_mass)
        self.ed_loads = QLineEdit()
        self.ed_loads.setPlaceholderText("sample µl per lane, e.g. 2:5, 3:5 (→ ng/µl)")
        f.addRow("sample loads", self.ed_loads)
        for w in (self.sp_ref, self.cb_mass, self.ed_loads):
            self._expert_items.append((f, w))
        b = QPushButton("Calibrate && quantify")
        b.clicked.connect(lambda: self.quantify())   # clicked passes "checked" as 1st arg
        f.addRow(b)
        self.lbl_cal = QLabel("")
        self.lbl_cal.setWordWrap(True)
        f.addRow(self.lbl_cal)
        lay.addWidget(g)
        self._preset_changed(DEFAULT_PRESET)

        # --- display
        g = QGroupBox("Compare lanes")
        f = QFormLayout(g)
        self.cb_norm = QComboBox()
        self.cb_norm.addItems(list(NORM_METHODS))
        self.cb_norm.setToolTip("Western blots: divide each band by the lane's loading-control band (e.g. actin) "
                                "or by the lane's total signal (total-protein / stain-free image).")
        f.addRow("normalise to", self.cb_norm)
        self.sp_lc = QDoubleSpinBox()
        self.sp_lc.setRange(0, 1e7)
        self.sp_lc.setDecimals(1)
        self.sp_lc.setSpecialValueText("–")
        self.sp_lc.setToolTip("Size of the loading-control band in the marker unit, e.g. 42 (kDa, actin)")
        f.addRow("loading-control size", self.sp_lc)
        self.sp_lc_tol = QDoubleSpinBox()
        self.sp_lc_tol.setRange(0.01, 0.5)
        self.sp_lc_tol.setSingleStep(0.01)
        self.sp_lc_tol.setValue(Params.lc_tol)
        f.addRow("  size tolerance", self.sp_lc_tol)
        self.cb_lc_ch = QComboBox()
        self.cb_lc_ch.setToolTip("Two-colour blots / separate total-protein image: where to read the loading control")
        f.addRow("  loading control from", self.cb_lc_ch)
        self.cb_tgt_ch = QComboBox()
        self.cb_tgt_ch.setToolTip("Two-colour blots: where to measure the bands that are normalised")
        f.addRow("  bands measured in", self.cb_tgt_ch)
        self.sp_ref_lane = QSpinBox()
        self.sp_ref_lane.setRange(0, 999)
        self.sp_ref_lane.setSpecialValueText("none")
        self.sp_ref_lane.setToolTip("Each band relative to the same band in this lane, e.g. the untreated "
                                    "control. Load the same reference sample on every gel to compare gels "
                                    "in a batch.")
        f.addRow("reference lane", self.sp_ref_lane)
        self.chk_rna = QCheckBox("RNA integrity (28S/18S, degradation)")
        self.chk_rna.setToolTip("Total RNA: 28S/18S ratio and the share of signal below 18S, per lane")
        f.addRow(self.chk_rna)
        for w in (self.sp_lc_tol, self.cb_lc_ch, self.cb_tgt_ch):
            self._expert_items.append((f, w))
        for w in (self.cb_norm, self.cb_lc_ch, self.cb_tgt_ch):
            w.currentIndexChanged.connect(self._compare_changed)
        for w in (self.sp_lc, self.sp_lc_tol, self.sp_ref_lane):
            w.valueChanged.connect(self._compare_changed)
        self.chk_rna.toggled.connect(self._compare_changed)
        lay.addWidget(g)

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
        self.grp_results = g
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
        b = QPushButton("Batch: more gels with these settings…")
        b.setToolTip("Analyse many images with this gel's settings (lanes, gel rows, ladder, marker "
                     "lanes, size range, detection). Each gel gets its own export and session; "
                     "batch_report.md lists the gels that need a look.")
        b.clicked.connect(lambda: self.batch())   # not batch(images=False)
        self.btn_batch = b
        f.addWidget(b)
        row = QHBoxLayout()
        b = QPushButton("Make sample sheet…")
        b.setToolTip("Gels with different layouts in one batch: choose the pictures, and autogel writes an "
                     "Excel sheet with one row per gel, filled in with the current settings. Change the "
                     "lanes, marker lanes, ladder or lane names where a gel differs.")
        b.clicked.connect(lambda: self.make_sample_sheet())
        self.btn_make_sheet = b
        row.addWidget(b)
        b = QPushButton("Batch with sample sheet…")
        b.setToolTip("Analyse the gels listed in a sample sheet, each with its own lanes, marker lanes, "
                     "ladder and lane names. Empty cells use the current settings.")
        b.clicked.connect(lambda: self.batch_from_sheet())
        self.btn_sheet_batch = b
        row.addWidget(b)
        f.addLayout(row)
        lay.addWidget(g)
        lay.addStretch()

        sc = QScrollArea()
        sc.setWidgetResizable(True)
        sc.setWidget(inner)
        self._scroll = sc
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(sc)
        # long ladder names must not widen the panel (the full name shows in the open list)
        for cb in inner.findChildren(QComboBox):
            cb.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            cb.setMinimumContentsLength(12)
            cb.view().setMinimumWidth(cb.view().sizeHintForColumn(0) + 30)
        # wide enough for the content in both modes: otherwise the dock opens narrower and cuts
        # off the right side (Help, Log, second buttons) behind a sideways scrollbar
        need = inner.minimumSizeHint().width() + sc.verticalScrollBar().sizeHint().width() + 4
        self.setMinimumWidth(max(380, need))

    # ------------------------------------------------------------ help
    def show_guide(self, section: str | None = None):
        if getattr(self, "_guide", None) is None:
            self._guide = HelpDialog(self)
        if section:
            self._guide.show_section(section)
        self._guide.show()
        self._guide.raise_()

    def show_tour(self):
        if getattr(self, "_tour", None) is not None:
            self._tour.close()
        self._tour = TourDialog(self, self)
        self._tour.show()
        self._tour.raise_()

    def shutdown(self):
        """Stop timers and detach from the analysis/viewer (closing the window or the dock)."""
        for t in (getattr(self, "_timer", None), getattr(self, "_zoom_timer", None)):
            if t is not None:
                t.stop()
        if self.ga is not None and self._on_log in self.ga.listeners:
            self.ga.listeners.remove(self._on_log)
        if getattr(self, "_zoom_hooked", False):
            with contextlib.suppress(Exception):
                self._camera().events.zoom.disconnect(self._zoom_cb)
            self._zoom_hooked = False

    def closeEvent(self, event):  # noqa: N802 (Qt API)
        self.shutdown()
        super().closeEvent(event)

    def _toggle_log(self, on: bool):
        self.txt_log.setVisible(on)
        self.btn_log.setText("Log ▾" if on else "Log ▸")

    def _on_log(self, e: dict):
        """Listener on GelAnalysis._log: every step appears in the status line and the log."""
        msg = describe(e)
        self.lbl_status.setText(f"✔ {msg}")
        self.txt_log.appendPlainText(f"{e['time'][11:]}  {msg}")
        self.viewer.status = msg

    @contextlib.contextmanager
    def _busy(self, what: str):
        """Show what is running (status line + wait cursor) for operations that take a moment."""
        from qtpy.QtWidgets import QApplication
        self.lbl_status.setText(f"⏳ {what} …")
        self.viewer.status = f"{what} …"
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()
        try:
            yield
        except Exception as e:
            self.lbl_status.setText(f"✖ {what} failed: {e}")
            raise
        finally:
            QApplication.restoreOverrideCursor()

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
        p.min_gel_rel = self.sp_gelrel.value()
        p.noise_floor = self.chk_floor.isChecked()
        p.deconvolve = self.chk_deconv.isChecked()
        p.split_shoulders = self.chk_shoulder.isChecked()
        p.checker_filter = self.sp_checker.value()
        p.n_lanes = self.sp_lanes.value() or None
        p.n_tiers = self.sp_tiers.value()
        p.n_columns = self.sp_cols.value()
        p.curved_lanes = self.chk_curved.isChecked()
        p.auto_crop = self.chk_crop.isChecked()
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
        p.n_tiers = self.sp_tiers.value()
        p.n_columns = self.sp_cols.value()
        if not self.ga:
            p.invert = None
        self.ga = None
        with self._busy(f"Loading {os.path.basename(path)}"):
            ga = GelAnalysis(path, p, orientation=orientation)
            self._set_analysis(ga)

    def _set_analysis(self, ga: GelAnalysis):
        self.ga = ga
        self.txt_log.clear()
        for e in ga.log:
            self.txt_log.appendPlainText(f"{e['time'][11:]}  {describe(e)}")
        if ga.log:
            self.lbl_status.setText(f"✔ {describe(ga.log[-1])}")
        ga.listeners.append(self._on_log)
        self.df = None
        self._focus_lane = None
        self._last_checked = None
        self._updating = True
        try:
            self.chk_invert.blockSignals(True)
            self.chk_invert.setChecked(bool(ga.gel.orientation.get("invert")))
            self.chk_invert.blockSignals(False)
            self.sp_angle.setValue(ga.gel.orientation.get("angle", 0.0))
            for name in ("gel", "saturated", "highlight", "analysed region", "gel splits", "lanes", "bands",
                         "sizes left", "sizes right"):
                if name in self.viewer.layers:
                    self.viewer.layers.remove(name)
            self._fill_channel_combos()
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
            region = self.viewer.add_shapes(name="analysed region", shape_type="rectangle",
                                            edge_color="#999999", face_color="transparent", edge_width=3)
            region.events.data.connect(lambda e: self._schedule("crop"))
            tiers = self.viewer.add_shapes(name="gel splits", shape_type="line", edge_color="#009E73",
                                           edge_width=4)
            tiers.events.data.connect(lambda e: self._schedule("tiers"))
            bands = self.viewer.add_shapes(name="bands", shape_type="rectangle",
                                           edge_color=self.colors["band"],
                                           face_color="transparent", edge_width=1.5)
            # size labels live on their own (non-editable) point layers beside the lanes
            for name, anchor in (("sizes left", "upper_right"), ("sizes right", "upper_left")):
                self.viewer.add_points(np.empty((0, 2)), name=name, size=1,
                                            face_color="transparent", border_color="transparent",
                                            border_width=0,
                                            text={"string": "{label}", "anchor": anchor,
                                                  "color": self.colors["label_text"], "size": 8},
                                            features={"label": np.empty(0, dtype=str)})
            self.viewer.layers.selection.active = bands
            lanes.events.data.connect(lambda e: self._schedule("lanes"))
            bands.events.data.connect(lambda e: self._schedule("bands"))
            lanes.events.highlight.connect(lambda e: self._lane_selected())
            bands.events.highlight.connect(lambda e: self._band_selected())
            self.viewer.reset_view()
            if not getattr(self, "_zoom_hooked", False):
                self._zoom_timer = QTimer(self)
                self._zoom_timer.setSingleShot(True)
                self._zoom_timer.setInterval(120)
                self._zoom_timer.timeout.connect(self._place_labels)
                self._zoom_cb = lambda e: self._zoom_timer.start()
                self._camera().events.zoom.connect(self._zoom_cb)
                self._zoom_hooked = True
        finally:
            self._updating = False
        exp = (ga.gel.info or {}).get("exposure")
        self.lbl_file.setText(os.path.basename(ga.path) + f"  ({a.shape[1]}×{a.shape[0]} px, {ga.gel.dtype}"
                              + (f", exposure {exp}" if exp else "") + ")")
        self.lbl_warn.setText("\n".join("⚠ " + w for w in ga.gel.warnings))
        self._refresh_layers()

    def _sync_controls_from_params(self):
        p = self.ga.params
        self.sp_lanes.setValue(p.n_lanes or 0)
        self.sp_tiers.blockSignals(True)
        self.sp_tiers.setValue(p.n_tiers or 1)
        self.sp_tiers.blockSignals(False)
        self.sp_cols.blockSignals(True)
        self.sp_cols.setValue(p.n_columns or 1)
        self.sp_cols.blockSignals(False)
        self.cb_mass.blockSignals(True)
        self.cb_mass.setCurrentText({v: k for k, v in MASS_MODELS.items()}.get(p.mass_model, "size-dependent"))
        self.cb_mass.blockSignals(False)
        self.cb_bg.setCurrentText(p.bg_method)
        self.sp_snr.setValue(p.min_snr)
        self.sp_cov.setValue(p.min_coverage)
        self.sp_gelrel.setValue(p.min_gel_rel)
        self.chk_floor.setChecked(p.noise_floor)
        for w, v in ((self.sp_lc, p.lc_size or 0), (self.sp_lc_tol, p.lc_tol),
                     (self.sp_ref_lane, p.reference_lane or 0)):
            w.blockSignals(True)
            w.setValue(v)
            w.blockSignals(False)
        self.cb_norm.blockSignals(True)
        self.cb_norm.setCurrentText({v: k for k, v in NORM_METHODS.items()}.get(p.norm_method, "none"))
        self.cb_norm.blockSignals(False)
        self.chk_rna.blockSignals(True)
        self.chk_rna.setChecked(p.rna_integrity)
        self.chk_rna.blockSignals(False)
        self._fill_channel_combos()
        self.chk_curved.setChecked(p.curved_lanes)
        self.chk_crop.setChecked(p.auto_crop or bool(p.crop))
        self.chk_deconv.setChecked(p.deconvolve)
        self.chk_shoulder.setChecked(p.split_shoulders)
        self.sp_checker.setValue(p.checker_filter)
        if p.marker_preset in PRESETS:
            self.cb_preset.setCurrentText(p.marker_preset)
        if p.marker_sizes:
            self.ed_sizes.setText(", ".join(f"{s:g}" for s in p.marker_sizes))
        self.ed_masses.setText(", ".join(f"{m:g}" for m in p.marker_masses) if p.marker_masses else "")
        self.ed_marker_lanes.setText(", ".join(map(str, p.marker_lanes)))
        for sp, v in ((self.sp_rmin, p.size_min or 0), (self.sp_rmax, p.size_max or 0)):
            sp.blockSignals(True)
            sp.setValue(v)
            sp.blockSignals(False)
        self.cb_range.blockSignals(True)
        self.cb_range.setCurrentText("custom" if (p.size_min or p.size_max) else "automatic (whole ladder)")
        self.cb_range.blockSignals(False)
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
        # the skew was measured in the old orientation: estimate it again (angle=None)
        self._reorient(lambda o: o.update(rot90_k=(o.get("rot90_k", 0) + 1) % 4, angle=None))

    def _flip_lr(self):
        # mirroring reverses the tilt
        self._reorient(lambda o: o.update(flip_lr=not o.get("flip_lr", False), angle=-(o.get("angle") or 0.0)))

    def _flip_ud(self):
        self._reorient(lambda o: o.update(flip_ud=not o.get("flip_ud", False), angle=-(o.get("angle") or 0.0)))

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
            with self._busy("Estimating lanes"):
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
        with self._busy(f"Detecting {p.n_lanes} lanes"):
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

    def _tiers_changed(self, n: int):
        if self.ga is None:
            return
        self.ga.params.n_tiers = n
        self.ga.params.tier_rows = None          # search the split(s) again
        self.ga._log("set_tiers", {"n_tiers": n, "reason": "gel rows setting"}, f"{n} gel row(s)")
        if self.ga.lanes:
            self.detect_lanes()

    def _cols_changed(self, n: int):
        if self.ga is None:
            return
        self.ga.params.n_columns = n
        self.ga.params.column_splits = None      # search the split(s) again
        self.ga._log("set_columns", {"n_columns": n, "reason": "gels side by side setting"},
                     f"{n} gel(s) side by side")
        if self.ga.lanes:
            self.detect_lanes()

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
        with self._busy("Detecting bands in all lanes"):
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
        lanes_changed = bands_changed = False
        if "lanes" in pending:
            layer = self.viewer.layers["lanes"]
            drawn = [(_lane_shape(l), l) for l in self.ga.lanes]
            lanes, unchanged = [], 0
            for s, kind in zip(layer.data, layer.shape_type):
                s = np.asarray(s, float)
                same = next((l for d, l in drawn if d.shape == s.shape and np.allclose(d, s, atol=0.5)), None)
                if same is not None:                 # untouched lane (incl. its bend)
                    lanes.append(same)
                    unchanged += 1
                    continue
                new = _lane_from_polygon(s) if kind == "polygon" else None
                if new is None:
                    y0, x0, y1, x1 = _bbox(s)
                    new = Lane(x=(x0 + x1) / 2, width=x1 - x0, y0=y0, y1=y1)
                # keep the name/source of the lane it replaces (same place at the top)
                prev = min(self.ga.lanes, key=lambda l: abs(l.x - new.x), default=None)
                if prev is not None and abs(prev.x - new.x) < prev.width / 2:
                    new.name, new.source = prev.name, "user"
                else:
                    new.source = "user"
                lanes.append(new)
            # napari can emit a data event without a real change: only real lane edits count
            if unchanged != len(self.ga.lanes) or len(lanes) != len(self.ga.lanes):
                self.ga.set_lanes(lanes, reason="edited in viewer")
                lanes_changed = True
        if "crop" in pending and "analysed region" in self.viewer.layers:
            shapes = self.viewer.layers["analysed region"].data
            new = None
            if len(shapes):
                y0, x0, y1, x1 = _bbox(shapes[0])
                new = [y0, y1, x0, x1]
            cur = list(self.ga.region())
            if new is not None and [round(v) for v in new] != [round(v) for v in cur]:
                self.ga.set_crop(new, reason="region moved in viewer")
                with self._busy("Re-detecting lanes and bands in the new region"):
                    self.ga.detect_lanes()
                    self.ga.detect_bands()
                lanes_changed = True
            elif new is None and self.ga.params.crop:
                self.ga.set_crop(None, reason="region deleted in viewer")
                lanes_changed = True
        if "tiers" in pending and "gel splits" in self.viewer.layers:
            rows, cols = [], []
            for ln in self.viewer.layers["gel splits"].data:
                ln = np.asarray(ln)
                dy, dx = np.ptp(ln[:, 0]), np.ptp(ln[:, 1])
                (rows if dx >= dy else cols).append(float(np.mean(ln[:, 0 if dx >= dy else 1])))
            changed = False
            p = self.ga.params
            if rows and [round(r) for r in sorted(rows)] != [round(r) for r in (p.tier_rows or [])]:
                self.ga.set_tier_rows(rows, reason="split line moved in viewer")
                changed = True
            if cols and [round(c) for c in sorted(cols)] != [round(c) for c in (p.column_splits or [])]:
                self.ga.set_column_splits(cols, reason="split line moved in viewer")
                changed = True
            if changed:
                for sp, v in ((self.sp_tiers, p.n_tiers), (self.sp_cols, p.n_columns)):
                    sp.blockSignals(True)
                    sp.setValue(v)
                    sp.blockSignals(False)
                with self._busy("Re-detecting lanes and bands for the new gel layout"):
                    self.ga.detect_lanes()
                    self.ga.detect_bands()
                lanes_changed = True
        if "bands" in pending and not lanes_changed:
            old = {(b.lane, round(b.y0), round(b.y1)): b for b in self.ga.bands}
            bands, orphan = [], 0
            for s in self.viewer.layers["bands"].data:
                y0, x0, y1, x1 = _bbox(s)
                y1 -= 1  # rectangles are drawn to b.y1 + 1 (inclusive last row)
                xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
                hit = [i for i, l in enumerate(self.ga.lanes) if l.contains(xc, yc)]
                if not hit:
                    orphan += 1
                    continue
                i = hit[0]
                prev = old.get((i, round(y0), round(y1)))
                bands.append(prev if prev else self.ga.band_for(i, y0, y1, source="user"))
            # napari also emits data events when the layer is redrawn or selected: only a real
            # change counts, otherwise redraw -> event -> "edit" -> quantify -> redraw ... loops
            if len(bands) != len(self.ga.bands) or any(b not in old.values() for b in bands):
                self.ga.set_bands(bands, reason="edited in viewer")
                bands_changed = True
            if orphan:
                self._info(f"{orphan} band rectangle(s) outside any lane were ignored")
        if not (lanes_changed or bands_changed):      # redraw events only: nothing to do
            return
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
                lanes.add([_lane_shape(l) for l in ga.lanes],
                          shape_type=["polygon" if l.bend else "rectangle" for l in ga.lanes],
                          edge_color=[col["marker"] if i in marker else col["lane"]
                                      for i in range(len(ga.lanes))])
                lanes.features = {"name": [l.name for l in ga.lanes]}
                # lane names centred ABOVE the lanes so they never sit on band labels
                h = float(np.median([l.y1 - l.y0 for l in ga.lanes]))
                lanes.text = {"string": "{name}", "anchor": "center", "color": col["lane_text"],
                              "size": 8, "translation": [-(h / 2 + 25), 0]}
            if "analysed region" in self.viewer.layers:
                rl = self.viewer.layers["analysed region"]
                rl.data = []
                H, W = ga.data.shape
                y0, y1, x0, x1 = ga.region()
                if (y0, y1, x0, x1) != (0, H, 0, W):
                    rl.add_rectangles([_rect(y0, x0, y1, x1)], edge_color="#999999", edge_width=3)
            if "gel splits" in self.viewer.layers:
                tl = self.viewer.layers["gel splits"]
                tl.data = []
                H, W = ga.data.shape
                rows = ga.params.tier_rows if (ga.params.n_tiers or 1) > 1 else []
                cols = ga.params.column_splits if (ga.params.n_columns or 1) > 1 else []
                lines = [np.array([[r, 0], [r, W - 1]], float) for r in (rows or [])]
                lines += [np.array([[0, c], [H - 1, c]], float) for c in (cols or [])]
                if lines:
                    tl.add_lines(lines, edge_color="#009E73", edge_width=4)
            bands = self.viewer.layers["bands"]
            bands.data = []
            if ga.bands:
                from matplotlib.colors import to_rgba
                edge = [col["band"] if b.source in ("auto", "ladder") else AI_BAND if b.source == "ai"
                        else col["user_band"] for b in ga.bands]
                # doubtful bands in grey: low band-checker probability or not confirmed by the AI
                doubt = set(self._doubtful())
                edge = [DOUBT_BAND if i in doubt else e for i, e in enumerate(edge)]
                # light fill: short (faint) bands stay visible as a box when zoomed out,
                # instead of their two edges merging into one line
                face = [to_rgba(c, 0.18) for c in edge]
                def box(b):
                    l, d = ga.lanes[b.lane], ga.lanes[b.lane].dx(b.y)
                    return _rect(b.y0, l.x0 + d, b.y1 + 1, l.x1 + d)
                bands.add_rectangles([box(b) for b in ga.bands], edge_color=edge, face_color=face)
        finally:
            self._updating = False
        self._update_band_labels()
        self._lane_selected()

    def _update_band_labels(self):
        """Size labels for marker lanes and the selected lane (or all bands when ticked).
        They sit beside the lane, never on the band; overlapping ones are hidden for the
        current zoom (reference bands first) and come back when you zoom in."""
        if self.ga is None or "sizes left" not in self.viewer.layers:
            return
        show = {L - 1 for L in self.ga.params.marker_lanes}
        if self._focus_lane is not None:
            show.add(self._focus_lane)
        if self.chk_all_labels.isChecked():
            show = set(range(len(self.ga.lanes)))
        self._labels = plan_labels(self.ga, self.df, show)
        self._place_labels()

    def _camera(self):
        scene = getattr(self.viewer, "scene", None)   # napari >= 0.9; older: viewer.camera
        return scene.camera if scene is not None and hasattr(scene, "camera") else self.viewer.camera

    def _place_labels(self):
        if self.ga is None or "sizes left" not in self.viewer.layers:
            return
        labels = getattr(self, "_labels", [])
        zoom = max(float(self._camera().zoom), 1e-6)         # screen px per image px
        font_px = 13.0                                         # 8 pt text incl. line spacing
        keep = declutter(labels, font_px / zoom)
        dy = 0.4 * font_px / zoom   # napari draws the text just above its anchor point: shift down
        for name, side in (("sizes left", "left"), ("sizes right", "right")):
            sel = [lab for lab, k in zip(labels, keep) if k and lab.side == side]
            layer = self.viewer.layers[name]
            self._updating = True
            try:
                layer.data = np.array([[lab.y + dy, lab.x] for lab in sel]).reshape(-1, 2)
                layer.features = {"label": [lab.text for lab in sel]}
                layer.text.color = self.colors["label_text"]
                layer.refresh_text()
            finally:
                self._updating = False

    def _check_keys(self) -> dict[int, tuple]:
        """Bands worth a look -> sort key, in the order 'Next band to check' visits them:
        grey bands autogel is unsure about (low band-checker probability or not confirmed by
        the AI second opinion; not in marker lanes), most doubtful first; then the AI's
        suggestions (green)."""
        ga, df = self.ga, self.df
        if ga is None:
            return {}
        keys = {i: (1, 0.0, b.lane, b.y) for i, b in enumerate(ga.bands) if b.source == "ai"}
        if df is None or len(df) != len(ga.bands) or not len(df):
            return keys
        doubt = np.zeros(len(df), bool)
        prob = np.ones(len(df))
        if "band_probability" in df:
            prob = np.nan_to_num(df["band_probability"].to_numpy(float), nan=1.0)
            doubt |= prob < DOUBT_PROBABILITY
        if "ai_support" in df:
            doubt |= ~df["ai_support"].astype(bool).to_numpy()
        markers = {L - 1 for L in ga.params.marker_lanes}   # checked by the calibration instead
        for i, b in enumerate(ga.bands):
            if doubt[i] and b.source == "auto" and b.lane not in markers:
                keys[i] = (0, float(prob[i]), b.lane, b.y)
        return keys

    def _doubtful(self) -> list[int]:
        """Grey (unsure) automatic bands, most doubtful first."""
        keys = self._check_keys()
        return sorted((i for i, k in keys.items() if k[0] == 0), key=keys.get)

    def next_band_to_check(self):
        """Select and centre the next grey or green band, so it can be kept, moved or deleted.
        Continues after the band shown last, so deleting it does not skip the next one."""
        if not self._need():
            return
        keys = self._check_keys()
        todo = sorted(keys, key=keys.get)
        if not todo:
            self._last_checked = None
            msg = "No bands left to check: no grey (unsure) or green (AI) bands."
            self.lbl_status.setText("✔ " + msg)
            self._info(msg)
            return
        last = getattr(self, "_last_checked", None)
        pos = next((n for n, i in enumerate(todo) if last is None or keys[i] > last), 0)
        i = todo[pos]
        self._last_checked = keys[i]
        b = self.ga.bands[i]
        l = self.ga.lanes[b.lane]
        bands = self.viewer.layers["bands"]
        self.viewer.layers.selection.active = bands
        bands.mode = "select"
        bands.selected_data = {i}
        try:   # centre on the band; zoom in (never out) to about five lanes across
            cam = self._camera()
            cam.center = (float(b.y), float((l.x0 + l.x1) / 2 + l.dx(b.y)))
            view = self.viewer.window._qt_viewer          # private: size of the image view
            lane_w = max(float(np.median([q.x1 - q.x0 for q in self.ga.lanes])), 1.0)
            target = min(view.width() / (5 * lane_w), view.height() / (self.ga.data.shape[0] / 3))
            if np.isfinite(target) and target > float(cam.zoom):
                cam.zoom = target
        except Exception:  # noqa: BLE001 - centring is a convenience only
            pass
        self._band_selected()
        kind = "unsure (grey)" if b.source == "auto" else "suggested by the AI (green)"
        msg = (f"Band {pos + 1} of {len(todo)} to check: lane {b.lane + 1}, {kind}. "
               "Keep it, move it, or press Delete to remove it.")
        self.lbl_status.setText("▸ " + msg)
        self._info(msg)

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
            o = l.outline()
            shapes.append(_rect(l.y0, o[:, 1].min() - 3, l.y1, o[:, 1].max() + 3))
            widths.append(3)
        if band is not None and 0 <= band < len(self.ga.bands):
            b = self.ga.bands[band]
            l = self.ga.lanes[b.lane]
            d = l.dx(b.y)
            shapes.append(_rect(b.y0 - 3, l.x0 + d - 6, b.y1 + 4, l.x1 + d + 6))
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

    def _range_preset_changed(self, name: str):
        if name in RESOLUTION_GUIDE:
            lo, hi, _unit = RESOLUTION_GUIDE[name]
        elif name.startswith("automatic"):
            lo, hi = 0, 0
        else:
            return            # custom: keep the values the user types
        for sp, v in ((self.sp_rmin, lo), (self.sp_rmax, hi)):
            sp.blockSignals(True)
            sp.setValue(v)
            sp.blockSignals(False)
        self._apply_range(reason=name)

    def _range_values_changed(self, _=None):
        self.cb_range.blockSignals(True)
        self.cb_range.setCurrentText("custom")
        self.cb_range.blockSignals(False)
        self._apply_range(reason="custom")

    def _apply_range(self, reason: str = ""):
        if self.ga is None:
            return
        self.ga.set_size_range(self.sp_rmin.value() or None, self.sp_rmax.value() or None, reason=reason)
        if self.ga.params.marker_lanes:
            self.quantify(silent=True)

    def _refresh_presets(self, select: str | None = None):
        cur = select or self.cb_preset.currentText()
        self.cb_preset.blockSignals(True)
        self.cb_preset.clear()
        self.cb_preset.addItems(list(PRESETS) + ["custom"])
        self.cb_preset.blockSignals(False)
        self.cb_preset.setCurrentText(cur if cur in PRESETS else "custom")
        self._preset_changed(self.cb_preset.currentText())

    def _open_ladders(self):
        from .ladder_dialog import LadderDialog
        dlg = LadderDialog(self, select=self.cb_preset.currentText())
        dlg.exec()
        self._refresh_presets(dlg.chosen)
        if dlg.chosen:
            self._info(f"ladder: {dlg.chosen}")

    def _suggest_ladder(self):
        """Rank library ladders against the bands of the (first) marker lane."""
        if not self._need():
            return
        txt = self.ed_marker_lanes.text().replace(";", ",").replace(" ", ",")
        try:
            lane = int(next(t for t in txt.split(",") if t.strip())) - 1
        except (StopIteration, ValueError):
            lane = 0
        if not 0 <= lane < len(self.ga.lanes):
            QMessageBox.information(self, "autogel", "Detect lanes first and enter the marker lane.")
            return
        vols = self.ga._volumes()
        bl = [(b.y, v) for b, v in zip(self.ga.bands, vols) if b.lane == lane]
        hits = suggest_ladders([y for y, _ in bl], [v for _, v in bl])
        if not hits:
            QMessageBox.information(self, "autogel", f"Lane {lane + 1} has {len(bl)} bands - too few to "
                                    "recognise a ladder. Detect bands first.")
            return
        labels = []
        for h in hits:
            also = len(h["names_same_sizes"]) - 1
            labels.append(f"{h['name']}   ({h['n_sizes']} sizes, {len(bl)} bands"
                          + (f"; same sizes as {also} other" if also else "") + ")")
        choice, ok = QInputDialog.getItem(self, "Suggest ladder",
                                          f"Best matches for the {len(bl)} bands in lane {lane + 1} "
                                          "(best first - check with your datasheet):", labels, 0, False)
        if ok and choice:
            name = hits[labels.index(choice)]["name"]
            self.cb_preset.setCurrentText(name)
            self.ga._log("suggest_ladder", {"lane": lane + 1, "ranking": [h["name"] for h in hits]},
                         f"chosen: {name}")

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
        self.ga.params.mass_model = MASS_MODELS[self.cb_mass.currentText()]
        self.ga.params.deconvolve = self.chk_deconv.isChecked()
        self._apply_compare()
        return True

    # ------------------------------------------------------------ channels / compare
    def _fill_channel_combos(self):
        """(Re)fill the channel choices from the analysis; keeps the current selections."""
        names = self.ga.channel_names() if self.ga else ["grey"]
        multi = self.ga is not None and len(self.ga.gel.channels) > 1
        p = self.ga.params if self.ga else Params()
        entries = [("default image", None)]
        if multi:
            entries += [(f"{k + 1}: {n}", k) for k, n in enumerate(names)] + [("all channels", -1)]
        same = [("same image", None)] + ([(f"{k + 1}: {n}", k) for k, n in enumerate(names)] if multi else [])
        for cb, items, cur in ((self.cb_channel, entries, p.channel), (self.cb_marker_ch, same, p.marker_channel),
                               (self.cb_lc_ch, same, p.lc_channel), (self.cb_tgt_ch, same, p.target_channel)):
            cb.blockSignals(True)
            cb.clear()
            for label, data in items:
                cb.addItem(label, data)
            idx = next((i for i in range(cb.count()) if cb.itemData(i) == cur), 0)
            cb.setCurrentIndex(idx)
            cb.setEnabled(cb.count() > 1)
            cb.blockSignals(False)

    def _show_image(self):
        if self.ga is None or "gel" not in self.viewer.layers:
            return
        a = self.ga.data
        lo, hi = np.percentile(a, [0.5, 99.7])
        layer = self.viewer.layers["gel"]
        layer.data = a
        layer.contrast_limits = (float(lo), float(max(hi, lo + 1e-6)))
        if "saturated" in self.viewer.layers:
            self.viewer.layers["saturated"].data = self.ga.gel.raw_max.astype(np.float32)

    def _channel_changed(self, _=None):
        if self.ga is None:
            return
        ch = self.cb_channel.currentData()
        if ch == self.ga.params.channel:
            return
        with self._busy("Switching channel"):
            self.ga.set_channel(ch)
            self._show_image()
        if self.ga.bands:
            self.quantify(silent=True)

    def _add_image_dialog(self):
        if not self._need():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Add image of the same gel / membrane", "",
                                              "Images (*.tif *.tiff *.png *.jpg *.jpeg *.bmp);;All (*)")
        if not path:
            return
        try:
            with self._busy(f"Adding {os.path.basename(path)}"):
                k = self.ga.add_image(path)
        except Exception as e:
            QMessageBox.warning(self, "autogel", str(e))
            return
        self._fill_channel_combos()
        self._info(f"added as channel {k + 1}. If it shows the marker, choose it under 'marker lanes from'.")

    def _marker_channel_changed(self, _=None):
        if self.ga is None:
            return
        ch = self.cb_marker_ch.currentData()
        if ch == self.ga.params.marker_channel:
            return
        with contextlib.suppress(ValueError):       # unreadable input is reported by quantify
            self._apply_marker()
        with self._busy("Reading the marker lanes from the chosen image"):
            self.ga.set_marker_channel(ch)
        self._refresh_layers()
        if self.ga.params.marker_lanes:
            self.quantify(silent=True)

    def ai_second_opinion(self):
        if not self._need():
            return
        if not self.ga.bands:
            QMessageBox.information(self, "autogel", "Detect bands first.")
            return
        from . import ai
        ok, why = ai.available()
        if not ok and "not installed" in why:
            QMessageBox.information(self, "autogel", "The AI second opinion is an optional extra.\n\n"
                                    "Install it with:  pip install autogel[ai]\n(onnxruntime, about 17 MB)")
            return
        if not ok:
            m = ai.MODEL
            if QMessageBox.question(
                    self, "AI model", f"Download GelGenie's band-recognition model once?\n\n"
                    f"{m['name']}, {m['size_mb']:.0f} MB, licence {m['licence']}\nfrom {m['source']}\n"
                    "It is stored in your .autogel folder and checked against its checksum.") != \
                    QMessageBox.StandardButton.Yes:
                return
            try:
                with self._busy("Downloading the AI model"):
                    ai.download_model()
            except Exception as e:
                QMessageBox.warning(self, "autogel", f"Download failed: {e}")
                return
        try:
            with self._busy("AI second opinion (GelGenie network)"):
                self.ga.ai_review()
        except Exception as e:
            QMessageBox.warning(self, "autogel", str(e))
            return
        self.quantify(silent=True)
        self._refresh_layers()

    def _apply_compare(self):
        if self.ga is None:
            return
        p = self.ga.params
        p.norm_method = NORM_METHODS[self.cb_norm.currentText()]
        p.lc_size = self.sp_lc.value() or None
        p.lc_tol = self.sp_lc_tol.value()
        p.lc_channel = self.cb_lc_ch.currentData()
        p.target_channel = self.cb_tgt_ch.currentData()
        p.reference_lane = self.sp_ref_lane.value() or None
        p.rna_integrity = self.chk_rna.isChecked()

    def _compare_changed(self, _=None):
        if self.ga is None or not self.ga.bands:
            return
        self.quantify(silent=True)

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
            how = "size-dependent" if mf["model"] == "size" else "through origin"
            msg.append(f"mass curve ({how}): typical error ±{mf['typical_error']:.0%} "
                       f"(leave-one-out), n={mf['n']} ladder bands")
        msg += getattr(self.ga, "compare_warnings", [])
        rna = getattr(self.ga, "rna", None)
        if rna is not None and len(rna):
            msg.append("RNA 28S/18S: " + ", ".join(
                f"{r.lane_name}: {r.rna_28s_18s:.2f} ({r.rna_below_18s_pct:.0f} % below 18S)"
                for r in rna.itertuples() if r.rna_28s_18s == r.rna_28s_18s))
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
            with self._busy("Exporting tables, profiles, image and report"):
                out = export_all(self.ga, d)
        except (ValueError, OSError) as e:   # bad input or file locked (e.g. xlsx open in Excel)
            QMessageBox.critical(self, "autogel", f"Export failed:\n{e}")
            return
        self._info(f"exported to {out}")

    def batch(self, images: list[str] | None = None, out_dir: str | None = None,
              sheet: list[dict] | None = None):
        if not self._need():
            return
        from .batch import run_batch
        try:
            self._apply_marker()
        except ValueError as e:
            QMessageBox.critical(self, "autogel", str(e))
            return
        if images is None and sheet is None:
            images, _ = QFileDialog.getOpenFileNames(
                self, "Gels to analyse with the current settings", os.path.dirname(self.ga.path),
                "Images (*.tif *.tiff *.png *.jpg *.jpeg *.bmp)")
        images = list(images or [])
        if not images and not sheet:
            return
        n = len(images) + len(sheet or [])
        if out_dir is None:
            first = images[0] if images else sheet[0]["path"]
            out_dir = QFileDialog.getExistingDirectory(self, "Output folder for the batch",
                                                       os.path.dirname(first))
        if not out_dir:
            return
        from qtpy.QtWidgets import QApplication, QProgressDialog
        dlg = QProgressDialog("Batch analysis …", "Cancel", 0, n, self)
        dlg.setWindowTitle("autogel batch")
        dlg.setMinimumDuration(0)

        def progress(i, n, name):
            dlg.setValue(i)
            dlg.setLabelText(f"Gel {i + 1} of {n}: {name}")
            self.lbl_status.setText(f"⏳ Batch: gel {i + 1} of {n} ({name}) …")
            QApplication.processEvents()
            return not dlg.wasCanceled()

        self.ga._log("batch", {"n_images": n, "out_dir": out_dir, "sample_sheet": bool(sheet)},
                     "started with these settings")
        summary = run_batch(images, self.ga.params, out_dir, progress, sheet=sheet)
        dlg.setValue(n)
        counts = summary["status"].value_counts().to_dict() if len(summary) else {}
        done = (f"Batch cancelled after {len(summary)} of {n} gel(s)" if summary.attrs.get("cancelled")
                else f"Batch finished: {len(summary)} gel(s)")
        msg = (f"{done} - {counts.get('ok', 0)} ok, "
               f"{counts.get('check', 0)} to check, {counts.get('error', 0)} error(s).\n"
               f"Overview: {os.path.join(out_dir, 'batch_report.md')}")
        self.lbl_status.setText("✔ " + msg.splitlines()[0])
        QMessageBox.information(self, "autogel batch", msg)
        return summary

    def make_sample_sheet(self, images: list[str] | None = None, path: str | None = None):
        """Write an Excel sheet (one row per chosen picture) filled in with the current settings."""
        if not self._need():
            return
        from .batch import write_sample_sheet
        try:
            self._apply_marker()
        except ValueError as e:
            QMessageBox.critical(self, "autogel", str(e))
            return
        if images is None:
            images, _ = QFileDialog.getOpenFileNames(
                self, "Pictures for the sample sheet", os.path.dirname(self.ga.path),
                "Images (*.tif *.tiff *.png *.jpg *.jpeg *.bmp)")
        if not images:
            return
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save sample sheet", os.path.join(os.path.dirname(images[0]), "sample_sheet.xlsx"),
                "Excel (*.xlsx)")
        if not path:
            return
        if not path.lower().endswith(".xlsx"):
            path += ".xlsx"
        try:
            write_sample_sheet(path, images, self.ga.params)
        except OSError as e:     # e.g. the sheet is open in Excel
            QMessageBox.critical(self, "autogel", f"Could not save the sample sheet:\n{e}")
            return
        msg = (f"Sample sheet saved: {path}\nOpen it in Excel, change the rows of gels that differ, save it, "
               "then click 'Batch with sample sheet…'.")
        self.lbl_status.setText(f"✔ Sample sheet saved ({len(images)} pictures)")
        QMessageBox.information(self, "autogel", msg)
        return path

    def batch_from_sheet(self, sheet_path: str | None = None, out_dir: str | None = None):
        """Batch over the gels of a sample sheet, each with its own lanes, ladder and names."""
        if not self._need():
            return
        from .batch import read_sample_sheet
        if sheet_path is None:
            sheet_path, _ = QFileDialog.getOpenFileName(
                self, "Sample sheet", os.path.dirname(self.ga.path), "Sample sheet (*.xlsx *.xls *.csv)")
        if not sheet_path:
            return
        try:
            sheet = read_sample_sheet(sheet_path)
        except (ValueError, OSError) as e:
            QMessageBox.critical(self, "autogel", str(e))
            return
        except Exception as e:  # noqa: BLE001 - unreadable file: say so plainly
            QMessageBox.critical(self, "autogel", f"Could not read the sample sheet ({type(e).__name__}): {e}")
            return
        return self.batch(images=[], out_dir=out_dir, sheet=sheet)

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

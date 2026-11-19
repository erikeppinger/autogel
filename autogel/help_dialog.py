"""In-app help: the plain-language user guide and a quick tour of the panel.

The guide (``help/guide.md``) ships with the program, so it works without internet. The full
technical manual and the issue form are on the project website and open in the browser.
"""
from __future__ import annotations

import contextlib
import os

from qtpy.QtCore import QUrl
from qtpy.QtGui import QDesktopServices
from qtpy.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from .markers import REPO

GUIDE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "help", "guide.md")
MANUAL_URL = f"https://github.com/{REPO}/blob/main/docs/USER_MANUAL.md"
REPORT_URL = f"https://github.com/{REPO}/issues/new"


def guide_sections(path: str = GUIDE_PATH) -> list[tuple[str, str]]:
    """[(title, markdown text)] of the guide, split at its top-level headings."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    out: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith("# "):
            out.append((line[2:].strip(), []))
        elif out:
            out[-1][1].append(line)
    return [(t, "\n".join(body).strip()) for t, body in out]


def open_url(url: str) -> bool:
    return bool(QDesktopServices.openUrl(QUrl(url)))


class HelpDialog(QDialog):
    """The user guide: sections on the left, text on the right, a search box on top."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("autogel – user guide")
        self.resize(820, 620)
        self.sections = guide_sections()
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search the guide (e.g. ladder, saturated, rotate)")
        self.search.textChanged.connect(self._filter)
        self.search.returnPressed.connect(self._find_next)
        row.addWidget(self.search, 1)
        b = QPushButton("Full manual (online)")
        b.setToolTip("Every setting in detail and how accurate autogel is (opens in your browser)")
        b.clicked.connect(lambda: open_url(MANUAL_URL))
        row.addWidget(b)
        lay.addLayout(row)
        row = QHBoxLayout()
        self.toc = QListWidget()
        self.toc.setMaximumWidth(220)
        self.toc.currentRowChanged.connect(self._show_row)
        row.addWidget(self.toc)
        self.text = QTextBrowser()
        self.text.setOpenExternalLinks(True)
        row.addWidget(self.text, 1)
        lay.addLayout(row, 1)
        self._visible = list(range(len(self.sections)))
        self._fill()

    def _fill(self):
        self.toc.blockSignals(True)
        self.toc.clear()
        for k in self._visible:
            self.toc.addItem(self.sections[k][0])
        self.toc.blockSignals(False)
        if self._visible:
            self.toc.setCurrentRow(0)
            self._show_row(0)
        else:
            self.text.setMarkdown("*Nothing found. Try another word.*")

    def _show_row(self, row: int):
        if 0 <= row < len(self._visible):
            title, body = self.sections[self._visible[row]]
            self.text.setMarkdown(f"# {title}\n\n{body}")
            q = self.search.text().strip()
            if q:
                self.text.find(q)

    def _filter(self, q: str):
        q = q.strip().lower()
        self._visible = [k for k, (t, b) in enumerate(self.sections) if not q or q in t.lower() or q in b.lower()]
        self._fill()

    def _find_next(self):
        q = self.search.text().strip()
        if q and not self.text.find(q):            # wrap to the next matching section
            row = self.toc.currentRow()
            if row + 1 < self.toc.count():
                self.toc.setCurrentRow(row + 1)
            elif self.toc.count():
                self.toc.setCurrentRow(0)

    def show_section(self, title: str):
        self.search.clear()
        for row, k in enumerate(self._visible):
            if self.sections[k][0] == title:
                self.toc.setCurrentRow(row)
                return


# (title, text, attribute of the autogel panel to point at or None)
TOUR = [
    ("Welcome",
     "This short tour shows you around the autogel panel. You work through it from top to "
     "bottom: open a picture, find the lanes, find the bands, choose your ladder, save the "
     "results.\n\nYou can leave the tour at any time and start it again from **Help**.", None),
    ("Easy or Expert",
     "**Easy** shows only what most gels need. **Expert** adds more settings. Your choice is "
     "remembered.", "btn_easy"),
    ("What is happening?",
     "This line always tells you what autogel is doing or has just done. **Log** shows every "
     "step for this gel – useful to retrace what you did.", "lbl_status"),
    ("1  Your picture",
     "Click **Open image…**. autogel turns the picture so the lanes run from top to bottom. "
     "If it is the wrong way round, use **Rotate 90°** or the flip buttons.\n\n"
     "Tip: use the original file from the imaging device (usually .tif), not a screenshot.",
     "grp_image"),
    ("2  Lanes",
     "Type the **number of lanes** (including ladder lanes) and click **Detect lanes**. "
     "Each lane gets a box.\n\nTo correct a lane, click **lanes** in the layer list on the left "
     "and move, resize, delete or draw boxes in the picture.", "grp_lanes"),
    ("3  Bands",
     "Click **Detect bands (all lanes)**. Grey boxes are bands autogel is unsure about. "
     "**AI second opinion** double-checks the bands on your computer.\n\nCorrect bands like "
     "lanes, in the layer **bands**.", "grp_bands"),
    ("4  Marker",
     "Choose your **ladder** and enter the **marker lane(s)**, then click "
     "**Calibrate & quantify**. autogel then knows the size of every band, and with the amounts "
     "from the ladder's datasheet roughly how much is in each band.\n\n"
     "Not sure which ladder? Try **Suggest ladder**.", "grp_marker"),
    ("5  Results",
     "**Export…** saves a table for Excel, the picture with all boxes, and a report of every step. "
     "**Save session…** keeps your work so you can continue later.", "grp_results"),
    ("Help is always here",
     "**Help** opens this tour again, the user guide (also works without internet), the full "
     "manual, and a form to report a problem or suggest an idea.", "btn_help"),
]


class TourDialog(QDialog):
    """Step-by-step tour; each step scrolls to and frames the part of the panel it explains."""

    FRAME = "border: 2px solid #E69F00; border-radius: 4px;"

    def __init__(self, panel: QWidget, parent=None):
        super().__init__(parent)
        self.panel = panel
        self.setWindowTitle("autogel – quick tour")
        self.setModal(False)
        self.resize(420, 260)
        self.step = 0
        self._framed: QWidget | None = None
        self._old_css = ""
        lay = QVBoxLayout(self)
        self.lbl_title = QLabel()
        self.lbl_title.setStyleSheet("font-size: 15px; font-weight: bold")
        lay.addWidget(self.lbl_title)
        self.lbl_text = QTextBrowser()
        self.lbl_text.setOpenExternalLinks(True)
        lay.addWidget(self.lbl_text, 1)
        row = QHBoxLayout()
        self.lbl_count = QLabel()
        row.addWidget(self.lbl_count)
        row.addStretch()
        self.btn_back = QPushButton("Back")
        self.btn_back.clicked.connect(lambda: self.go(self.step - 1))
        self.btn_next = QPushButton("Next")
        self.btn_next.clicked.connect(self._next)
        row.addWidget(self.btn_back)
        row.addWidget(self.btn_next)
        lay.addLayout(row)
        self.go(0)

    def _next(self):
        if self.step + 1 < len(TOUR):
            self.go(self.step + 1)
        else:
            self.close()

    def go(self, k: int):
        k = max(0, min(len(TOUR) - 1, k))
        self.step = k
        title, text, target = TOUR[k]
        self.lbl_title.setText(title)
        self.lbl_text.setMarkdown(text)
        self.lbl_count.setText(f"{k + 1} of {len(TOUR)}")
        self.btn_back.setEnabled(k > 0)
        self.btn_next.setText("Finish" if k == len(TOUR) - 1 else "Next")
        self._frame(getattr(self.panel, target, None) if target else None)

    def _frame(self, w: QWidget | None):
        if self._framed is not None:
            with contextlib.suppress(RuntimeError):  # widget already deleted
                self._framed.setStyleSheet(self._old_css)
        self._framed = w
        if w is None:
            return
        self._old_css = w.styleSheet()
        name = w.objectName() or f"tour_{id(w)}"
        w.setObjectName(name)
        w.setStyleSheet(self._old_css + f" #{name} {{ {self.FRAME} }}")
        scroll = getattr(self.panel, "_scroll", None)
        if scroll is not None:
            scroll.ensureWidgetVisible(w, 20, 20)

    def closeEvent(self, e):
        self._frame(None)
        super().closeEvent(e)

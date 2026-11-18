"""Ladder manager: view, add, edit, delete, import/export and share ladders.

Built-in and ApE ladders are read-only; editing one and saving creates *your* copy (user
library, ``~/.autogel/ladders.json``). Sharing opens a pre-filled GitHub issue in the browser;
nothing is sent without the user clicking.
"""
from __future__ import annotations

from qtpy.QtCore import Qt, QUrl
from qtpy.QtGui import QDesktopServices
from qtpy.QtWidgets import (
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from . import markers

ORIGIN = {"builtin": "built-in, checked against vendor data", "update": "library update from GitHub",
          "ape": "imported from ApE (sizes only, not checked)", "user": "your ladder"}


def _numbers(text: str) -> list[float]:
    return [float(t) for t in text.replace(";", ",").replace(" ", ",").split(",") if t.strip()]


class LadderDialog(QDialog):
    def __init__(self, parent=None, select: str | None = None):
        super().__init__(parent)
        self.setWindowTitle("Ladders")
        self.resize(820, 520)
        self.chosen: str | None = None          # ladder to use when the dialog closes
        root = QHBoxLayout(self)
        left = QVBoxLayout()
        self.ed_filter = QLineEdit()
        self.ed_filter.setPlaceholderText("search (name, vendor, catalog)")
        self.ed_filter.textChanged.connect(self._fill)
        left.addWidget(self.ed_filter)
        self.lst = QListWidget()
        self.lst.currentTextChanged.connect(self._show)
        left.addWidget(self.lst)
        root.addLayout(left, 2)

        right = QVBoxLayout()
        self.lbl_origin = QLabel("")
        self.lbl_origin.setWordWrap(True)
        right.addWidget(self.lbl_origin)
        f = QFormLayout()
        self.ed_name, self.ed_vendor, self.ed_catalog = QLineEdit(), QLineEdit(), QLineEdit()
        self.cb_unit = QComboBox()
        self.cb_unit.addItems(list(markers.UNITS))
        self.ed_sizes = QLineEdit()
        self.ed_sizes.setPlaceholderText("e.g. 10000, 8000, 6000 … (any order, pasted from the datasheet)")
        self.ed_masses = QLineEdit()
        self.ed_masses.setPlaceholderText("optional: ng (DNA) / µg (protein) per band, same order as sizes")
        self.sp_ref = QDoubleSpinBox()
        self.sp_ref.setRange(0, 1e6)
        self.sp_ref.setDecimals(3)
        self.sp_ref.setSpecialValueText("–")
        self.cb_ref_unit = QComboBox()
        self.cb_ref_unit.addItems(["µg", "µl", "ng"])
        ref = QHBoxLayout()
        ref.addWidget(self.sp_ref)
        ref.addWidget(self.cb_ref_unit)
        self.ed_source = QLineEdit()
        self.ed_source.setPlaceholderText("https://… datasheet")
        self.ed_note = QLineEdit()
        for label, w in (("name", self.ed_name), ("vendor", self.ed_vendor), ("catalog no.", self.ed_catalog),
                         ("unit", self.cb_unit), ("sizes", self.ed_sizes), ("masses", self.ed_masses)):
            f.addRow(label, w)
        f.addRow("masses refer to", ref)
        f.addRow("source", self.ed_source)
        f.addRow("note", self.ed_note)
        right.addLayout(f)
        self.lbl_msg = QLabel("")
        self.lbl_msg.setWordWrap(True)
        right.addWidget(self.lbl_msg)

        def row(*buttons):
            h = QHBoxLayout()
            for txt, fn, tip in buttons:
                b = QPushButton(txt)
                b.setToolTip(tip)
                b.clicked.connect(fn)
                h.addWidget(b)
            right.addLayout(h)
        row(("New", self._new, "Empty form for a new ladder"),
            ("Save as my ladder", self._save, "Save to your library (a copy, if this is a built-in ladder)"),
            ("Delete", self._delete, "Delete one of your own ladders"))
        row(("Import…", self._import, "Add ladders from a file (e.g. from a colleague)"),
            ("Export…", self._export, "Save this ladder to a file to share it"),
            ("Suggest for library…", self._share,
             "Open a pre-filled GitHub issue so it can be checked and shipped with the next release"))
        row(("Update library from GitHub…", self._update,
             "Download the newest curated ladder list (between releases)"),
            ("Use this ladder", self._use, "Select it in the panel and close"))
        right.addStretch(1)
        root.addLayout(right, 3)
        self._fill()
        if select and select in markers.PRESETS:
            items = self.lst.findItems(select, Qt.MatchFlag.MatchExactly)
            if items:
                self.lst.setCurrentItem(items[0])

    # ------------------------------------------------------------ list / form
    def _fill(self, *_):
        q = self.ed_filter.text().lower()
        cur = self.lst.currentItem().text() if self.lst.currentItem() else None
        self.lst.blockSignals(True)
        self.lst.clear()
        for n, e in markers.PRESETS.items():
            hay = " ".join(str(e.get(k, "")) for k in ("name", "vendor", "catalog")).lower()
            if q in hay:
                self.lst.addItem(n)
        self.lst.blockSignals(False)
        if cur:
            items = self.lst.findItems(cur, Qt.MatchFlag.MatchExactly)
            if items:
                self.lst.setCurrentItem(items[0])

    def _show(self, name: str):
        e = markers.PRESETS.get(name)
        if not e:
            return
        self.ed_name.setText(e["name"])
        self.ed_vendor.setText(e.get("vendor") or "")
        self.ed_catalog.setText(e.get("catalog") or "")
        self.cb_unit.setCurrentText(e.get("unit", "bp"))
        self.ed_sizes.setText(", ".join(f"{s:g}" for s in e["sizes"]))
        self.ed_masses.setText(", ".join(f"{m:g}" for m in e["masses"]) if e.get("masses") else "")
        ref = e.get("reference_load") or {}
        self.sp_ref.setValue(float(ref.get("amount") or 0))
        self.cb_ref_unit.setCurrentText(ref.get("unit") or "µg")
        self.ed_source.setText(e.get("source") or "")
        self.ed_note.setText(e.get("note") or "")
        origin = ORIGIN.get(e.get("library", ""), "")
        extra = "" if e.get("user_defined") else " - edit and 'Save as my ladder' to make your own copy"
        self.lbl_origin.setText(f"<b>{name}</b><br>{origin}{extra}")
        self.lbl_msg.setText("")

    def _entry(self) -> dict:
        sizes = _numbers(self.ed_sizes.text())
        masses = _numbers(self.ed_masses.text()) or None
        ref = ({"amount": self.sp_ref.value(), "unit": self.cb_ref_unit.currentText(),
                "text": f"masses for {self.sp_ref.value():g} {self.cb_ref_unit.currentText()}"}
               if self.sp_ref.value() > 0 else None)
        return {"name": self.ed_name.text().strip(), "vendor": self.ed_vendor.text().strip(),
                "catalog": self.ed_catalog.text().strip(), "unit": self.cb_unit.currentText(),
                "sizes": sizes, "masses": masses, "reference_load": ref,
                "source": self.ed_source.text().strip(), "note": self.ed_note.text().strip()}

    def _say(self, msg: str, ok: bool = True):
        self.lbl_msg.setStyleSheet("" if ok else "color: orange")
        self.lbl_msg.setText(msg)

    # ---------------------------------------------------------------- actions
    def _new(self):
        self.lst.clearSelection()
        for w in (self.ed_name, self.ed_vendor, self.ed_catalog, self.ed_sizes, self.ed_masses,
                  self.ed_source, self.ed_note):
            w.clear()
        self.sp_ref.setValue(0)
        self.lbl_origin.setText("<b>new ladder</b>")
        self.ed_name.setFocus()

    def _save(self) -> bool:
        try:
            e = self._entry()
            existing = markers.PRESETS.get(e["name"])
            if existing and not existing.get("user_defined"):
                e["name"] += " (my copy)"
            markers.save_user_ladder(e["name"], e["sizes"], e["unit"], e["masses"], source=e["source"],
                                     note=e["note"], vendor=e["vendor"], catalog=e["catalog"],
                                     reference_load=e["reference_load"])
        except ValueError as err:
            self._say(f"Not saved: {err}", ok=False)
            return False
        self._fill()
        items = self.lst.findItems(e["name"], Qt.MatchFlag.MatchExactly)
        if items:
            self.lst.setCurrentItem(items[0])
        self._say(f"Saved '{e['name']}' to your library.")
        return True

    def _delete(self):
        name = self.ed_name.text().strip()
        e = markers.PRESETS.get(name)
        if not e or not e.get("user_defined"):
            self._say("Only your own ladders can be deleted.", ok=False)
            return
        if QMessageBox.question(self, "Delete ladder", f"Delete your ladder '{name}'?") != \
                QMessageBox.StandardButton.Yes:
            return
        markers.delete_user_ladder(name)
        self._fill()
        self._new()
        self._say(f"Deleted '{name}'.")

    def _import(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import ladders", "", "Ladder files (*.json);;All (*)")
        if not path:
            return
        r = markers.import_ladders(path)
        if r["skipped"] and QMessageBox.question(
                self, "Import ladders", "Already in the library: " + ", ".join(r["skipped"])
                + "\nReplace them with the imported version?") == QMessageBox.StandardButton.Yes:
            r2 = markers.import_ladders(path, overwrite=True)
            r["replaced"], r["skipped"] = r2["replaced"], []
        self._fill()
        parts = [f"{len(r['added'])} added", f"{len(r['replaced'])} replaced"]
        if r["skipped"]:
            parts.append(f"{len(r['skipped'])} skipped")
        msg = "Import: " + ", ".join(parts)
        if r["errors"]:
            msg += ". Problems: " + " | ".join(r["errors"][:3])
        self._say(msg, ok=not r["errors"])

    def _export(self):
        name = self.ed_name.text().strip()
        if name not in markers.PRESETS:
            self._say("Save the ladder first.", ok=False)
            return
        safe = "".join(c if c.isalnum() or c in "-_ " else "_" for c in name).strip()
        path, _ = QFileDialog.getSaveFileName(self, "Export ladder", f"{safe}.json", "Ladder files (*.json)")
        if path:
            markers.export_ladders([name], path)
            self._say(f"Exported to {path}")

    def _share(self):
        e = self._entry()
        errs = markers.validate_ladder(e)
        if errs:
            self._say("Fix first: " + "; ".join(errs), ok=False)
            return
        if not e["source"]:
            self._say("Add the datasheet link (source) so the entry can be checked.", ok=False)
            return
        if QMessageBox.question(
                self, "Suggest for library",
                "This opens a pre-filled GitHub issue in your browser (a free GitHub account is needed "
                "to submit it). Nothing is sent until you submit it there.\n\nOpen it now?") != \
                QMessageBox.StandardButton.Yes:
            return
        QDesktopServices.openUrl(QUrl(markers.contribution_url(e)))
        self._say("Opened in your browser. Once checked, it ships with the next release and "
                  "reaches everyone via 'Update library from GitHub'.")

    def _update(self):
        if QMessageBox.question(
                self, "Update library",
                f"Download the current ladder list from github.com/{markers.REPO}?\n"
                "Your own ladders are kept.") != QMessageBox.StandardButton.Yes:
            return
        try:
            n, _ = markers.update_library()
        except Exception as err:            # network, proxy, invalid file - keep the old library
            self._say(f"Update failed, old library kept: {err}", ok=False)
            return
        self._fill()
        self._say(f"Library updated: {n} curated ladders.")

    def _use(self):
        """Use the ladder shown in the form; unsaved edits are saved first (as your ladder)."""
        try:
            e = self._entry()
        except ValueError as err:
            self._say(f"Check the numbers: {err}", ok=False)
            return
        stored = markers.PRESETS.get(e["name"])
        same = (stored is not None
                and [float(x) for x in stored["sizes"]] == sorted(e["sizes"], reverse=True)
                and [float(x) for x in (stored.get("masses") or [])] == [float(x) for x in (e["masses"] or [])])
        if not same:
            if not self._save():
                return
            e["name"] = self.lst.currentItem().text()
        self.chosen = e["name"]
        self.accept()

"""Consistency checks on the built-in ladder library."""
import json
import os

import pytest

import autogel

with open(os.path.join(os.path.dirname(autogel.__file__), "ladders.json"), encoding="utf-8") as _f:
    LIB = json.load(_f)
LADDERS = LIB["ladders"]


@pytest.mark.parametrize("e", LADDERS, ids=[e["name"] for e in LADDERS])
def test_entry_is_consistent(e):
    sizes = e["sizes"]
    assert e["unit"] in ("bp", "nt", "kDa")
    assert sizes == sorted(sizes, reverse=True) and len(set(sizes)) == len(sizes)
    if e.get("masses"):
        assert len(e["masses"]) == len(sizes)
    for r in e.get("reference_bands", []):
        assert r in sizes
    if e.get("verified"):
        assert e.get("source", "").startswith("https://") and e.get("checked")


@pytest.mark.parametrize("e", [e for e in LADDERS if e.get("masses") and e.get("reference_load")],
                         ids=lambda e: e["name"])
def test_masses_add_up_to_reference_load(e):
    """Datasheet band masses must add up to the stated load (catches transcription errors)."""
    total_ng = sum(e["masses"])
    ref_ng = e["reference_load"]["amount"] * 1000
    assert total_ng == pytest.approx(ref_ng, rel=0.02)


def test_default_preset_exists():
    from autogel.markers import DEFAULT_PRESET, PRESETS
    assert DEFAULT_PRESET in PRESETS and PRESETS[DEFAULT_PRESET]["verified"]


def test_ape_import_agrees_with_curated_entries():
    """ApE's NEB lists (doublets merged) must equal the sizes checked against NEB's datasheets."""
    path = os.path.join(os.path.dirname(autogel.__file__), "ladders_ape.json")
    with open(path, encoding="utf-8") as f:
        ape = {e["catalog"]: e for e in json.load(f)["ladders"]}
    curated = {e["catalog"]: e for e in LADDERS if e.get("verified") and e.get("catalog")}
    common = set(ape) & set(curated)
    assert {"N3200", "N3231", "N3232"} <= common
    for c in common:
        assert ape[c]["sizes"] == curated[c]["sizes"], c


def test_library_merges_ape_without_duplicates():
    from autogel.markers import load_library
    lib = load_library()
    cats = [e["catalog"] for e in lib.values() if e.get("catalog")]
    assert len(cats) == len(set(cats))
    assert len(lib) > 60

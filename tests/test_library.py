"""User ladder library: validation, import/export, delete, update, sharing link, suggestions."""
import json
from urllib.parse import parse_qs, urlparse

import numpy as np
import pytest

from autogel import markers


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("AUTOGEL_HOME", str(tmp_path / "home"))
    markers._reload()
    yield tmp_path
    monkeypatch.delenv("AUTOGEL_HOME")
    markers._reload()


GOOD = {"name": "Lab ladder", "unit": "bp", "sizes": [3000, 1000, 500, 200], "masses": [30, 20, 10, 5],
        "reference_load": {"amount": 0.065, "unit": "µg"}, "source": "https://example.org/ds.pdf"}


@pytest.mark.parametrize("bad, msg", [
    ({**GOOD, "unit": "kb"}, "unit"),
    ({**GOOD, "sizes": [100, 1000]}, "top to bottom"),
    ({**GOOD, "sizes": [1000, 1000, 5]}, "top to bottom"),
    ({**GOOD, "sizes": [1000]}, "at least 2"),
    ({**GOOD, "sizes": [1000, -5]}, "positive"),
    ({**GOOD, "sizes": [1000, "x"]}, "positive"),
    ({**GOOD, "sizes": [1000, float("nan")]}, "positive"),
    ({**GOOD, "masses": [1, 2]}, "masses"),
    ({**GOOD, "masses": [1, -2, 3, 4]}, "non-negative"),
    ({**GOOD, "name": " "}, "name"),
    ({**GOOD, "reference_load": {"amount": 0}}, "reference_load"),
    ({**GOOD, "source": 5}, "source"),
    ("not a dict", "object"),
])
def test_validate_rejects(bad, msg):
    errs = markers.validate_ladder(bad)
    assert errs and any(msg in e for e in errs)


def test_validate_accepts_library():
    assert markers.validate_ladder(GOOD) == []
    for e in markers.load_library().values():
        assert markers.validate_ladder(e) == [], e["name"]


def test_save_sorts_pairs_and_delete(home):
    markers.save_user_ladder("mine", [500, 3000, 1000], "bp", [10, 30, 20])
    e = markers.PRESETS["mine"]
    assert e["sizes"] == [3000, 1000, 500] and e["masses"] == [30, 20, 10]     # masses stay with sizes
    assert e["user_defined"] and not e["verified"]
    assert markers.delete_user_ladder("mine") and "mine" not in markers.PRESETS
    assert not markers.delete_user_ladder(markers.DEFAULT_PRESET)               # built-in: refused
    assert markers.DEFAULT_PRESET in markers.PRESETS


def test_export_import_roundtrip(home):
    markers.save_user_ladder("mine", [3000, 1000, 500], "bp", [30, 20, 10], source="https://x.org")
    f = home / "share.json"
    assert markers.export_ladders(["mine", "does not exist"], str(f)) == 1
    data = json.loads(f.read_text(encoding="utf-8"))
    assert "user_defined" not in data["ladders"][0]
    markers.delete_user_ladder("mine")
    r = markers.import_ladders(str(f))
    assert r["added"] == ["mine"] and markers.PRESETS["mine"]["masses"] == [30, 20, 10]
    r = markers.import_ladders(str(f))                       # again: exists -> skipped
    assert r["skipped"] == ["mine"] and not r["added"]
    r = markers.import_ladders(str(f), overwrite=True)
    assert r["replaced"] == ["mine"]


def test_import_single_object_bad_entries_and_garbage(home):
    one = home / "one.json"
    one.write_text(json.dumps(GOOD), encoding="utf-8")
    assert markers.import_ladders(str(one))["added"] == ["Lab ladder"]
    mixed = home / "mixed.json"
    mixed.write_text(json.dumps({"ladders": [{**GOOD, "name": "ok2"}, {**GOOD, "name": "bad", "unit": "x"}]}),
                     encoding="utf-8")
    r = markers.import_ladders(str(mixed))
    assert r["added"] == ["ok2"] and len(r["errors"]) == 1 and "bad" in r["errors"][0]
    junk = home / "junk.json"
    junk.write_text("{not json", encoding="utf-8")
    assert markers.import_ladders(str(junk))["errors"]
    assert markers.import_ladders(str(home / "missing.json"))["errors"]
    empty = home / "empty.json"
    empty.write_text("[]", encoding="utf-8")
    assert markers.import_ladders(str(empty))["errors"]


def test_import_cannot_shadow_builtin_without_consent(home):
    f = home / "x.json"
    f.write_text(json.dumps({**GOOD, "name": markers.DEFAULT_PRESET}), encoding="utf-8")
    assert markers.import_ladders(str(f))["skipped"] == [markers.DEFAULT_PRESET]
    assert markers.PRESETS[markers.DEFAULT_PRESET]["verified"]


def test_broken_user_file_does_not_break_library(home):
    (home / "home").mkdir()
    (home / "home" / "ladders.json").write_text("{broken", encoding="utf-8")
    assert markers.DEFAULT_PRESET in markers.load_library()
    (home / "home" / "ladders.json").write_text(json.dumps({"ladders": [{"name": "bad"}, GOOD]}),
                                                encoding="utf-8")
    lib = markers.load_library()
    assert "Lab ladder" in lib and "bad" not in lib


def test_library_update_validates_and_overrides(home):
    with open(markers.__file__.replace("markers.py", "ladders.json"), encoding="utf-8") as f:
        lib = json.load(f)
    lib["ladders"].append({**GOOD, "name": "New vendor ladder", "verified": True})
    n, path = markers._store_update(json.dumps(lib).encode(), "https://example.org/l.json", "now")
    assert n == len(lib["ladders"]) and markers.PRESETS["New vendor ladder"]["library"] == "update"
    lib["ladders"].append({**GOOD, "name": "broken", "sizes": [1, 2]})
    with pytest.raises(ValueError, match="invalid"):
        markers._store_update(json.dumps(lib).encode(), "u", "now")
    assert "New vendor ladder" in markers.PRESETS                     # old update kept
    with pytest.raises(ValueError):
        markers._store_update(b"x" * (markers.MAX_LIBRARY_BYTES + 1), "u", "now")
    with pytest.raises(ValueError):
        markers.update_library("http://insecure.example.org/l.json")


def test_contribution_url_prefills_issue_form():
    url = markers.contribution_url({**GOOD, "user_defined": True, "library": "user"})
    u = urlparse(url)
    assert u.netloc == "github.com" and u.path.endswith("/issues/new")
    q = {k: v[0] for k, v in parse_qs(u.query).items()}
    assert q["template"] == "ladder.yml" and q["sizes"] == "3000, 1000, 500, 200"
    assert q["reference_load"] == "0.065 µg"
    assert "user_defined" not in json.loads(q["json"])


def _lane(sizes, intensity, rng, missing_last=False, curve=20.0):
    ls = np.log10(sizes)
    x = (ls - ls.min()) / (ls.max() - ls.min())
    ys = 150 + (1 - x) * 900 + curve * (1 - x) ** 2 + rng.normal(0, 1.5, len(x))
    vols = np.asarray(intensity, float) * rng.uniform(0.85, 1.15, len(x))
    n = len(x) - 1 if missing_last else len(x)
    return list(ys[:n]), list(vols[:n])


@pytest.mark.parametrize("name", [markers.DEFAULT_PRESET, "GeneRuler 1 kb Plus (SM1331)",
                                  "PageRuler Plus Prestained (26619)"])
def test_suggest_finds_ladder(name):
    if name not in markers.PRESETS:
        pytest.skip("ladder not in library")
    e = markers.PRESETS[name]
    rng = np.random.default_rng(1)
    for missing in (False, True):
        ys, vols = _lane(e["sizes"], e.get("masses") or [1] * len(e["sizes"]), rng, missing)
        hits = markers.suggest_ladders(ys, vols)
        assert any(name in h["names_same_sizes"] for h in hits), (name, missing, hits)


def test_suggest_edge_cases():
    assert markers.suggest_ladders([100, 200, 300]) == []              # too few bands
    hits = markers.suggest_ladders([100, 200, 300, 400, 500, 600], [1] * 6, unit="kDa")
    assert all(markers.PRESETS[h["name"]]["unit"] == "kDa" for h in hits)
    assert len(hits) <= 3


def test_cli_import_export(home, capsys):
    from autogel.cli import main
    f = home / "one.json"
    f.write_text(json.dumps(GOOD), encoding="utf-8")
    main(["ladders", "--import", str(f)])
    assert "added 1" in capsys.readouterr().out
    main(["ladders", "--import", str(f)])
    assert "use --replace" in capsys.readouterr().out
    out = home / "out.json"
    main(["ladders", "--export", "Lab ladder", str(out)])
    assert json.loads(out.read_text(encoding="utf-8"))["ladders"][0]["sizes"] == GOOD["sizes"]
    with pytest.raises(SystemExit):
        main(["ladders", "--export", "nope", str(out)])

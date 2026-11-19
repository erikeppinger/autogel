"""Band checker (learned) and AI second opinion (GelGenie network via onnxruntime)."""
import numpy as np
import pytest

from autogel.analysis import GelAnalysis, Params
from autogel.checker import FEATURES, band_features, load_model, predict

from .synth import ladder, make_gel, y_for_size

SIZES = [10000, 6000, 4000, 3000, 2000, 1500, 1000, 700, 500, 300, 100]


def _gel(tmp_path, **kw):
    lad, _ = ladder(sizes=SIZES)
    lanes = {0: lad, 9: lad}
    for k in range(1, 9):
        lanes[k] = [(y_for_size(2500), 6000, 5), (y_for_size(800), 3000, 5)]
    return make_gel(str(tmp_path / "g.tif"), lanes=lanes, noise=40, **kw)


def test_checker_model_ships_and_scores(tmp_path):
    m = load_model()
    assert m is not None and m["features"] == FEATURES and len(m["coef"]) == len(FEATURES)
    assert 0.4 < m["test_precision_all"] < m["test_precision_at_0.2"]           # learned something
    ga = GelAnalysis(_gel(tmp_path), Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    f = band_features(ga)
    assert len(f) == len(ga.bands) and list(f.columns) == FEATURES
    p = predict(m, f)
    assert ((p >= 0) & (p <= 1)).all() and np.median(p) > 0.5                  # clean gel: mostly real
    df = ga.quantify()
    assert "band_probability" in df


def test_checker_filter_never_touches_marker_lanes(tmp_path):
    ga = GelAnalysis(_gel(tmp_path), Params(n_lanes=10, marker_lanes=[1, 10], checker_filter=0.95))
    ga.detect_lanes()
    ga.detect_bands()
    assert sum(b.lane == 0 for b in ga.bands) == len(SIZES)
    assert sum(b.lane == 9 for b in ga.bands) == len(SIZES)
    assert "unlikely (band checker)" in ga.log[-1]["result"]
    assert sum(b.lane not in (0, 9) for b in ga.bands) < 16                      # some samples removed


def _fake_network(monkeypatch, ga, keep_lanes):
    """Stand-in for the network: band pixels exactly where the synthetic bands of keep_lanes are,
    plus one band autogel cannot see (lane 5, at 5 kb)."""
    from autogel import ai
    m = np.zeros(ga.data.shape, bool)
    for b in ga.bands:
        if b.lane in keep_lanes:
            L = ga.lanes[b.lane]
            m[int(b.y0):int(b.y1) + 1, int(L.x0):int(L.x1)] = True
    L = ga.lanes[5]
    y = int(y_for_size(5000))
    m[y - 4:y + 5, int(L.x0) + 3:int(L.x1) - 3] = True
    monkeypatch.setattr(ai, "available", lambda: (True, "ok"))
    monkeypatch.setattr(ai, "oriented_mask", lambda ga_, threshold=0.5: m)


def test_ai_review_flags_and_adds(tmp_path, monkeypatch):
    ga = GelAnalysis(_gel(tmp_path), Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    n0 = len(ga.bands)
    _fake_network(monkeypatch, ga, keep_lanes={0, 1, 2, 3, 4, 5, 6, 7, 9})       # lane 9 (index 8) unconfirmed
    r = ga.ai_review()
    assert len(ga.bands) == n0 + 1 and sum(b.source == "ai" for b in ga.bands) == 1
    added = next(b for b in ga.bands if b.source == "ai")
    assert added.lane == 5 and abs(added.y - y_for_size(5000)) < 6
    un = [b for b, s in zip(ga.bands, r["support"]) if not s]
    assert un and all(b.lane == 8 for b in un)
    df = ga.quantify()
    assert "ai_support" in df and (~df[df.lane == 9]["ai_support"]).all()
    assert ga.log[-1]["step"] in ("quantify", "ai_review")
    # after an edit the stored opinion no longer matches the bands -> no stale column
    ga.set_bands(ga.bands[:-1])
    assert "ai_support" not in ga.quantify()


def test_ai_unavailable_is_explained(tmp_path, monkeypatch):
    from autogel import ai
    monkeypatch.setattr(ai, "available", lambda: (False, "the AI extra is not installed"))
    ga = GelAnalysis(_gel(tmp_path), Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    with pytest.raises(RuntimeError, match="not installed"):
        ga.ai_review()


def test_real_network_if_installed(tmp_path):
    from autogel import ai
    ok, why = ai.available()
    if not ok:
        pytest.skip(why)
    ga = GelAnalysis(_gel(tmp_path), Params(n_lanes=10))
    ga.detect_lanes()
    ga.detect_bands()
    r = ga.ai_review(add_missing=False)
    assert sum(r["support"]) >= 0.8 * len(ga.bands)                                # clean synthetic gel

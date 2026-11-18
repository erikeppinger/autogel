"""Regression tests on the lab's example gels (skipped when the share is not mounted)."""
import os

import pytest

from autogel.analysis import GelAnalysis, Params
from autogel.markers import DEFAULT_PRESET as PRESET

GEL_DIR = os.environ.get("AUTOGEL_TEST_GELS", "")


@pytest.mark.parametrize("name", ["Test_1.tif", "Test_2.tif"])
def test_example_gel(name):
    path = os.path.join(GEL_DIR, name)
    if not os.path.exists(path):
        pytest.skip("example gels not available")
    ga = GelAnalysis(path, Params(n_lanes=20))
    ga.detect_lanes()
    ga.detect_bands()
    assert len(ga.lanes) == 20
    assert sum(b.lane == 0 for b in ga.bands) == 19     # NEB 1 kb Plus: 19 bands
    assert sum(b.lane == 19 for b in ga.bands) == 19
    ga.set_marker([1, 20], preset=PRESET)
    df = ga.quantify()
    assert not any(c.warnings for c in ga.calibrations)
    assert df["size_bp"].notna().all()
    # mass curve: size-dependent response; a ladder band left out is predicted within ~25 %
    # (one slope through the origin: ~80 %, see lab book 0.4.0)
    assert ga.mass_fit["model"] == "size" and ga.mass_fit["typical_error"] < 0.3
    # ladder recognition: the 1 kb Plus pattern (bright 3 kb / 1 kb / 500 bp bands) ranks first
    from autogel.markers import suggest_ladders
    vols = ga._volumes()
    bl = [(b.y, v) for b, v in zip(ga.bands, vols) if b.lane == 0]
    hits = suggest_ladders([y for y, _ in bl], [v for _, v in bl])
    assert PRESET in hits[0]["names_same_sizes"], hits

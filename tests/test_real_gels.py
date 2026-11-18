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

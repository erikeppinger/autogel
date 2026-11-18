"""Regression test over the simulated gel gallery (tests/gallery.py).

Thresholds document what autogel achieves today. Known limitations are stated explicitly
instead of being hidden: smiling gels need a middle marker, JPEG/scratched gels give a few
spurious extra bands, saturated bands are under-quantified (and flagged).
"""
import pytest

from .gallery import SCENARIOS, evaluate

BY_NAME = {s.name: s for s in SCENARIOS}

# scenario -> (max median size error %, max extra bands, max median ratio error %)
EXPECT = {
    "agarose_clean": (1, 0, 5),
    "agarose_low_signal": (1, 0, 5),
    "agarose_uneven_light": (1, 0, 5),
    "agarose_tilted": (1, 0, 5),
    "agarose_smiling": (12, 0, 5),             # limitation: 2 edge markers cannot model a smile
    "agarose_smiling_3_markers": (4, 0, 5),    # ...a middle marker fixes most of it
    "agarose_two_combs": (1, 0, 5),            # two stacked gels, split found automatically
    "agarose_side_by_side": (1, 0, 5),         # two gels next to each other, split found automatically
    "agarose_dirty": (1, 4, 5),                # scratches crossing lanes can look like bands
    "agarose_smears_empty": (1, 0, 5),
    "agarose_overexposed": (2, 0, 25),         # clipped bands under-quantified (flagged saturated)
    "agarose_jpeg_8bit": (1, 14, 5),           # JPEG/8-bit artefacts -> spurious faint bands
    "agarose_rotated": (1, 0, 5),
    "agarose_many_narrow": (1, 0, 5),
    "agarose_few_wide": (1, 0, 5),
    "sds_coomassie": (1, 0, 5),
    "sds_coomassie_grey": (1, 0, 5),
}


def test_every_scenario_has_expectations():
    assert set(EXPECT) == set(BY_NAME)


@pytest.mark.parametrize("name", list(EXPECT))
def test_scenario(name, tmp_path):
    sc = BY_NAME[name]
    r = evaluate(sc, str(tmp_path))
    size_max, extra_max, ratio_max = EXPECT[name]
    assert r["lanes_found"] == r["lanes_expected"]
    assert r["wells_on_top"]
    assert r["bands_missed"] == 0
    assert r["median_size_error_pct"] <= size_max
    assert r["extra_bands"] <= extra_max
    if r["median_ratio_error_pct"] is not None:
        assert r["median_ratio_error_pct"] <= ratio_max


def test_overexposed_bands_are_flagged(tmp_path):
    from autogel.analysis import GelAnalysis, Params

    from .gallery import render
    sc = BY_NAME["agarose_overexposed"]
    path, _ = render(sc, str(tmp_path))
    ga = GelAnalysis(path, Params(n_lanes=sc.n_lanes))
    ga.detect_lanes()
    ga.detect_bands()
    df = ga.quantify()
    assert df.saturated.any() and any("saturated" in w for w in ga.gel.warnings)

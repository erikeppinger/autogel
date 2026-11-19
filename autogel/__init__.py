"""autogel - automatic lane/band detection and quantification for agarose and SDS gels."""
__version__ = "0.5.0"

from .analysis import Band, GelAnalysis, Lane, Params  # noqa: E402,F401
from .markers import PRESETS  # noqa: E402,F401

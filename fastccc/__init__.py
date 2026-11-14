import sys
import importlib.metadata
from loguru import logger

__version__ = importlib.metadata.version("fastccc")


LOG_FORMAT = '<cyan>{time:YYYY-MM-DD HH:mm:ss}</cyan> | <level>{level: <8}</level> | <level>{message}</level>'
_log_handler_id = None


def _configure_logging():
    """Install FastCCC's log format unless the caller already configured loguru.

    Only loguru's built-in default handler (id 0) is replaced; handlers added by
    the user are left untouched so library users keep control of logging.
    """
    global _log_handler_id
    try:
        logger.remove(0)
    except ValueError:          # default handler already removed → user has their own setup
        return
    if not logger._core.handlers:
        _log_handler_id = logger.add(sys.stdout, level="INFO", format=LOG_FORMAT)


def log_to(stream) -> None:
    """Send FastCCC's own log messages to ``stream`` (the command line uses stderr so that
    ``--json`` output on stdout stays machine-readable). Does nothing when the caller
    configured loguru themselves."""
    global _log_handler_id
    if _log_handler_id is None:
        return
    try:
        logger.remove(_log_handler_id)
    except ValueError:
        pass
    _log_handler_id = logger.add(stream, level="INFO", format=LOG_FORMAT)


_configure_logging()

# ── original entry points (unchanged) ────────────────────────────────────────
from .core import Cauchy_combination_of_statistical_analysis_methods
from .core import statistical_analysis_method

from . import build_reference
from . import infer_query
from .report import (
    generate_report,
    generate_infer_report,
    generate_reference_report,
)

# ── one-call API, inspection, resources and result objects ───────────────────
from .api import run, build_reference_panel, compare_to_reference, report_split, suggest_celltype_mapping
from .celltype_match import MappingSuggestion
from .datacheck import inspect_data as inspect, detect_matrix_type, suggest_celltype_key
from .resources import (
    list_databases, resolve_database,
    list_reference_panels, download_reference_panel, resolve_reference, fastccc_home,
)
from .results import (
    RunResult, MultiRunResult, ReferenceResult, ReferencePanel, load_result,
)

__all__ = [
    '__version__',
    'run', 'build_reference_panel', 'compare_to_reference', 'report_split', 'suggest_celltype_mapping', 'MappingSuggestion',
    'inspect', 'detect_matrix_type', 'suggest_celltype_key',
    'list_databases', 'resolve_database', 'list_reference_panels', 'download_reference_panel',
    'resolve_reference', 'fastccc_home',
    'RunResult', 'MultiRunResult', 'ReferenceResult', 'ReferencePanel', 'load_result',
    'Cauchy_combination_of_statistical_analysis_methods', 'statistical_analysis_method',
    'build_reference', 'infer_query',
    'generate_report', 'generate_infer_report', 'generate_reference_report',
]

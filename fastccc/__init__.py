import sys
from loguru import logger

logger.remove()
logger.add(sys.stdout, level="INFO", format='<cyan>{time:YYYY-MM-DD HH:mm:ss}</cyan> | <level>{level: <8}</level> | <level>{message}</level>')

from .core import Cauchy_combination_of_statistical_analysis_methods
from .core import statistical_analysis_method

from . import build_reference
from . import infer_query
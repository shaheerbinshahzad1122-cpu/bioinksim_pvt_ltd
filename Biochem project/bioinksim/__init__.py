"""
BioInkSim: A Biochemical and Rheological Simulator for Polymer-Based Tissue Engineering and 3D Bioprinting
===========================================================================================================
"""

from .data import (
    POLYMERS_DB,
    CELLS_DB,
    GROWTH_FACTORS_DB,
    NOZZLES_DB,
    CROSSLINKERS_DB,
    TISSUES_DB,
    BENCHMARKS_DB,
)
from .engine import BioInkEngine
from .calibration import RheologyCalibrator
from .uq import UncertaintyEngine
from .sensitivity import SensitivityEngine
from .report_pdf import generate_pdf_report

__version__ = "2.5.0"
__all__ = [
    "POLYMERS_DB",
    "CELLS_DB",
    "GROWTH_FACTORS_DB",
    "NOZZLES_DB",
    "CROSSLINKERS_DB",
    "TISSUES_DB",
    "BENCHMARKS_DB",
    "BioInkEngine",
    "RheologyCalibrator",
    "UncertaintyEngine",
    "SensitivityEngine",
    "generate_pdf_report",
]

"""
Tests for Polymer Rheology Calculations
=======================================
"""

from bioinksim.engine import BioInkEngine


def test_alginate_shear_thinning():
    """Alginate must be shear-thinning (n < 1.0) and have positive consistency index."""
    rheo = BioInkEngine.calculate_rheology("Alginate", 3.0, 180, 25.0, 7.4)
    assert rheo["n"] < 1.0
    assert rheo["K"] > 0
    assert "Shear-Thinning" in rheo["classification"]


def test_gelatin_temperature_transition():
    """
    Gelatin below gelation temp (22°C < 27.5°C) should have much higher consistency
    and yield stress than in the sol state (37°C).
    """
    rheo_gel = BioInkEngine.calculate_rheology("Gelatin", 5.0, 100, 22.0, 7.4)
    rheo_sol = BioInkEngine.calculate_rheology("Gelatin", 5.0, 100, 37.0, 7.4)

    assert rheo_gel["K"] > rheo_sol["K"] * 3.0
    assert rheo_gel["yield_stress_pa"] > rheo_sol["yield_stress_pa"] * 5.0


def test_chitosan_ph_precipitation():
    """Chitosan at pH > 6.5 forms aggregates leading to elevated apparent yield stress."""
    rheo_acid = BioInkEngine.calculate_rheology("Chitosan", 2.5, 150, 25.0, 5.5)
    rheo_base = BioInkEngine.calculate_rheology("Chitosan", 2.5, 150, 25.0, 7.5)

    assert rheo_base["yield_stress_pa"] > rheo_acid["yield_stress_pa"]

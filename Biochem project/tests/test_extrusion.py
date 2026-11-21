"""
Tests for Extrusion Mechanics, Bagley Correction, and Yield Stress Cutoff
=========================================================================
"""

import math
from bioinksim.engine import BioInkEngine
from bioinksim.data import NOZZLES_DB


def test_bagley_correction_ratio():
    """
    Verify that Bagley entrance loss correction does NOT divide by 2.4.
    For a 21G blunt needle (L/D ~ 24.7), delta_P_cap must be ~90-95% of delta_P_total,
    not 42%!
    """
    pressure_kpa = 150.0
    ext = BioInkEngine.calculate_extrusion_mechanics(
        nozzle_gauge="21G",
        nozzle_type="Cylindrical (Blunt)",
        pressure_kpa=pressure_kpa,
        syringe_vol_ml=5.0,
        K=10.0,
        n=0.5,
        yield_stress_pa=0.0,
    )

    delta_P_total = pressure_kpa * 1000.0
    delta_P_cap = ext["effective_capillary_pressure_pa"]
    ratio = delta_P_cap / delta_P_total

    # Ratio must be between 0.90 and 0.96
    assert 0.90 <= ratio <= 0.96, f"Bagley ratio for 21G needle should be ~0.93, got {ratio:.3f}"

    # Wall shear stress for 21G (R = 0.257 mm, L = 12.7 mm):
    # tau_w = (delta_P_cap * R) / (2 * L) ~ (139800 * 0.257e-3) / (2 * 0.0127) ~ 1414 Pa
    assert 1350 <= ext["wall_shear_stress_pa"] <= 1450, (
        f"Wall shear stress should be ~1414 Pa, got {ext['wall_shear_stress_pa']} Pa"
    )


def test_conical_vs_blunt_entrance_loss():
    """Conical needle should have lower entrance loss than blunt needle."""
    ext_blunt = BioInkEngine.calculate_extrusion_mechanics(
        "21G", "Cylindrical (Blunt)", 100.0, 5.0, 10.0, 0.5, 0.0
    )
    ext_conical = BioInkEngine.calculate_extrusion_mechanics(
        "21G", "Conical (Tapered)", 100.0, 5.0, 10.0, 0.5, 0.0
    )

    assert ext_conical["entrance_loss_pa"] < ext_blunt["entrance_loss_pa"], (
        "Conical tip must produce lower entrance pressure loss than blunt tip."
    )


def test_yield_stress_no_flow():
    """
    When tau_w <= yield_stress, no flow must occur (Q = 0, v = 0, flow_active = False).
    """
    # At 2.0 kPa on 21G, tau_w is ~ 18.8 Pa.
    # Set yield_stress = 30.0 Pa.
    ext = BioInkEngine.calculate_extrusion_mechanics(
        nozzle_gauge="21G",
        nozzle_type="Cylindrical (Blunt)",
        pressure_kpa=2.0,
        syringe_vol_ml=5.0,
        K=10.0,
        n=0.5,
        yield_stress_pa=30.0,
    )

    assert ext["flow_active"] is False, "Flow should be inactive when tau_w <= yield_stress"
    assert ext["flow_rate_mm3_s"] == 0.0, f"Flow rate must be 0.0, got {ext['flow_rate_mm3_s']}"
    assert ext["extrusion_velocity_mm_s"] == 0.0, f"Velocity must be 0.0, got {ext['extrusion_velocity_mm_s']}"
    assert "NO FLOW" in ext["flow_status"], f"Expected NO FLOW status, got {ext['flow_status']}"


def test_yield_stress_flow_when_exceeded():
    """When tau_w > yield_stress, material flows with Herschel-Bulkley mechanics."""
    ext = BioInkEngine.calculate_extrusion_mechanics(
        nozzle_gauge="21G",
        nozzle_type="Cylindrical (Blunt)",
        pressure_kpa=100.0,
        syringe_vol_ml=5.0,
        K=10.0,
        n=0.5,
        yield_stress_pa=30.0,
    )

    assert ext["flow_active"] is True
    assert ext["flow_rate_mm3_s"] > 0.0
    assert ext["extrusion_velocity_mm_s"] > 0.0


def test_syringe_volume_metrics():
    """Verify that syringe_vol_ml is utilized for plunger velocity and dispensing duration."""
    ext = BioInkEngine.calculate_extrusion_mechanics(
        nozzle_gauge="21G",
        nozzle_type="Cylindrical (Blunt)",
        pressure_kpa=150.0,
        syringe_vol_ml=5.0,
        K=10.0,
        n=0.5,
        yield_stress_pa=0.0,
    )

    assert "plunger_speed_mm_s" in ext
    assert ext["plunger_speed_mm_s"] > 0.0
    assert "dispensing_time_s" in ext
    assert ext["dispensing_time_s"] > 0.0

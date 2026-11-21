"""
Closed-Form Sanity Checks and Analytical Limit Verifications
===========================================================
"""

import math
from bioinksim.engine import BioInkEngine
from bioinksim.data import NOZZLES_DB


def test_newtonian_limit_reproduces_hagen_poiseuille():
    """
    In the Newtonian limit (n = 1.0, yield_stress = 0.0), the Rabinowitsch-Mooney
    capillary flow solution must EXACTLY reproduce the closed-form Hagen-Poiseuille equation:
    Q = (pi * R^4 * Delta_P_cap) / (8 * mu * L)
    """
    mu = 0.85  # Dynamic viscosity in Pa·s
    pressure_kpa = 120.0
    nozzle = "21G"
    nozzle_type = "Cylindrical (Blunt)"
    L_eff = 0.0127  # m
    D_i = NOZZLES_DB[nozzle]["inner_diam_mm"] * 1e-3
    R = D_i / 2.0

    ext = BioInkEngine.calculate_extrusion_mechanics(
        nozzle_gauge=nozzle,
        nozzle_type=nozzle_type,
        pressure_kpa=pressure_kpa,
        syringe_vol_ml=5.0,
        K=mu,
        n=1.0,
        yield_stress_pa=0.0,
    )

    delta_P_cap = ext["effective_capillary_pressure_pa"]

    # Exact Hagen-Poiseuille volumetric flow rate:
    Q_poiseuille_m3_s = (math.pi * (R ** 4) * delta_P_cap) / (8.0 * mu * L_eff)
    Q_poiseuille_mm3_s = Q_poiseuille_m3_s * 1e9

    # Check relative discrepancy
    rel_error = abs(ext["flow_rate_mm3_s"] - Q_poiseuille_mm3_s) / Q_poiseuille_mm3_s
    assert rel_error < 0.001, (
        f"Newtonian flow rate {ext['flow_rate_mm3_s']:.4f} diverges from "
        f"Hagen-Poiseuille {Q_poiseuille_mm3_s:.4f} by {rel_error*100:.3f}%"
    )


def test_zero_crosslinking_sol_state():
    """At zero crosslinker, the network cannot form and E must be ~0."""
    cross = BioInkEngine.calculate_crosslinking(
        crosslinker_name="CaCl2",
        conc_mm=0.0,
        time_min=0.0,
        temp_c=25.0,
        polymer_name="Alginate",
        conc_w_v=3.0,
    )

    assert cross["crosslinking_degree_pct"] == 0.0
    assert cross["youngs_modulus_kpa"] <= 0.1
    assert "Sol State" in cross["gel_state"]


def test_viability_bounded():
    """Viability must always remain physically bounded between 0 and 100%."""
    # Extreme shear stress condition (e.g. 30G needle at 500 kPa)
    ext = BioInkEngine.calculate_extrusion_mechanics("30G", "Cylindrical (Blunt)", 500.0, 5.0, 50.0, 0.4, 0.0)
    safety = BioInkEngine.calculate_cell_safety("Stem Cells", 20.0, ext["wall_shear_stress_pa"], ext["residence_time_s"], 4.5)

    assert 0.0 <= safety["viability_pct"] <= 100.0

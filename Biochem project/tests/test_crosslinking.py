"""
Tests for Hydrogel Crosslinking, Concentration Scaling, and Percolation Threshold
==================================================================================
"""

from bioinksim.engine import BioInkEngine


def test_modulus_scales_with_concentration():
    """
    Young's modulus must scale with polymer concentration as E ~ c^1.7.
    A 6% GelMA hydrogel must have significantly higher modulus than a 3% GelMA hydrogel.
    """
    cross_3pct = BioInkEngine.calculate_crosslinking(
        crosslinker_name="LAP + 405nm",
        conc_mm=1.7,
        time_min=1.0,
        temp_c=25.0,
        polymer_name="GelMA",
        conc_w_v=3.0,
    )
    cross_6pct = BioInkEngine.calculate_crosslinking(
        crosslinker_name="LAP + 405nm",
        conc_mm=1.7,
        time_min=1.0,
        temp_c=25.0,
        polymer_name="GelMA",
        conc_w_v=6.0,
    )

    ratio = cross_6pct["youngs_modulus_kpa"] / cross_3pct["youngs_modulus_kpa"]
    expected_ratio = (6.0 / 3.0) ** 1.7  # 2^1.7 ~ 3.25

    assert 3.0 <= ratio <= 3.5, (
        f"Expected 2-fold concentration increase to yield ~3.25x modulus, got {ratio:.2f}x"
    )


def test_percolation_threshold_gel_point():
    """
    Below the critical crosslinking percolation threshold chi_c, the hydrogel is in
    the sol state with G' ~ 0 and E ~ 0.
    Above chi_c, elasticity emerges sharply.
    """
    # chi_c for Alginate is 0.12 (12%)
    # Test very short crosslinking time (0.01 min) with low conc (0.5 mM) -> chi < chi_c
    cross_sol = BioInkEngine.calculate_crosslinking(
        crosslinker_name="CaCl2",
        conc_mm=0.5,
        time_min=0.05,
        temp_c=25.0,
        polymer_name="Alginate",
        conc_w_v=3.0,
    )

    assert cross_sol["crosslinking_degree_pct"] <= 12.0
    assert cross_sol["youngs_modulus_kpa"] <= 0.5, (
        f"Modulus below gel point should be near 0, got {cross_sol['youngs_modulus_kpa']} kPa"
    )
    assert "Sol State" in cross_sol["gel_state"]

    # Now crosslink with standard 50 mM for 10 min -> chi well above chi_c
    cross_gel = BioInkEngine.calculate_crosslinking(
        crosslinker_name="CaCl2",
        conc_mm=50.0,
        time_min=10.0,
        temp_c=25.0,
        polymer_name="Alginate",
        conc_w_v=3.0,
    )

    assert cross_gel["crosslinking_degree_pct"] > 50.0
    assert cross_gel["youngs_modulus_kpa"] > 50.0, (
        f"Modulus above gel point must be significant, got {cross_gel['youngs_modulus_kpa']} kPa"
    )


def test_ph_impact_on_alginate_crosslinking():
    """At acidic pH (pH < 4.0), alginate carboxylate protonation inhibits Ca2+ crosslinking."""
    cross_ph7 = BioInkEngine.calculate_crosslinking("CaCl2", 50.0, 10.0, 25.0, "Alginate", 3.0, ph=7.4)
    cross_ph3 = BioInkEngine.calculate_crosslinking("CaCl2", 50.0, 10.0, 25.0, "Alginate", 3.0, ph=3.0)

    assert cross_ph3["crosslinking_degree_pct"] < cross_ph7["crosslinking_degree_pct"]
    assert cross_ph3["youngs_modulus_kpa"] < cross_ph7["youngs_modulus_kpa"]

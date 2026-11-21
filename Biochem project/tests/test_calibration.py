"""
Tests for Rheology Calibration Curve Fitting
============================================
"""

from bioinksim.calibration import RheologyCalibrator


def test_power_law_fit_accuracy():
    """Verify that power law regression correctly recovers K and n."""
    true_K = 5.4
    true_n = 0.52
    rates = [0.5, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0, 500.0]
    points = [(g, true_K * (g ** true_n)) for g in rates]

    fit = RheologyCalibrator.fit_power_law(points)

    assert abs(fit["K"] - true_K) < 0.05, f"Expected K={true_K}, got {fit['K']}"
    assert abs(fit["n"] - true_n) < 0.01, f"Expected n={true_n}, got {fit['n']}"
    assert fit["r_squared"] > 0.999


def test_herschel_bulkley_fit():
    """Verify Herschel-Bulkley fitting recovers yield stress, K, and n."""
    true_tau_y = 12.0
    true_K = 7.5
    true_n = 0.45
    rates = [0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0]
    points = [(g, true_tau_y + true_K * (g ** true_n)) for g in rates]

    fit = RheologyCalibrator.fit_herschel_bulkley(points)

    assert abs(fit["yield_stress_pa"] - true_tau_y) < 1.0, f"Expected tau_y ~ {true_tau_y}, got {fit['yield_stress_pa']}"
    assert abs(fit["n"] - true_n) < 0.05, f"Expected n ~ {true_n}, got {fit['n']}"
    assert fit["r_squared"] > 0.99

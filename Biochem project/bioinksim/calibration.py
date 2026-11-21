"""
BioInkSim Rheology Calibration & Curve Fitting Studio
=====================================================
Fits experimental shear sweep data (shear rate vs shear stress / viscosity)
to non-Newtonian flow models:
1. Ostwald-de Waele Power-Law Model: tau = K * gamma_dot^n
2. Herschel-Bulkley Yield-Stress Model: tau = tau_y + K * gamma_dot^n

Calculates goodness-of-fit metrics (R^2, RMSE), parameter standard errors,
and generates parity curves for user inspection.
"""

import math
import csv
import io


class RheologyCalibrator:
    """Pure-Python non-linear regression engine for bioink rheometry."""

    @staticmethod
    def parse_data(raw_text, data_type="stress"):
        """
        Parses CSV or whitespace-delimited text containing experimental rheometer points.
        Format per line: shear_rate, shear_stress (or viscosity)
        """
        points = []
        # Strip comments and empty lines
        lines = [line.strip() for line in raw_text.strip().splitlines() if line.strip() and not line.strip().startswith("#")]
        
        # Check for header
        first_line = lines[0].lower()
        if any(h in first_line for h in ["rate", "gamma", "shear", "stress", "visc"]):
            lines = lines[1:]

        for line in lines:
            # Handle comma, tab, or space separation
            delimiter = "," if "," in line else ("\t" if "\t" in line else None)
            parts = [p.strip() for p in line.split(delimiter) if p.strip()]
            if len(parts) >= 2:
                try:
                    gamma_dot = float(parts[0])
                    val = float(parts[1])
                    if gamma_dot > 0 and val > 0:
                        if data_type == "viscosity":
                            # val is viscosity eta (Pa·s) -> tau = eta * gamma_dot
                            stress = val * gamma_dot
                        else:
                            stress = val
                        points.append((gamma_dot, stress))
                except ValueError:
                    continue

        points.sort(key=lambda pt: pt[0])
        return points

    @staticmethod
    def fit_power_law(points):
        """
        Fits Ostwald-de Waele model: tau = K * gamma_dot^n
        via linear regression in log-log space: ln(tau) = ln(K) + n * ln(gamma_dot).
        """
        if len(points) < 3:
            raise ValueError("At least 3 valid experimental data points are required for curve fitting.")

        n_pts = len(points)
        x = [math.log(pt[0]) for pt in points]
        y = [math.log(pt[1]) for pt in points]

        x_mean = sum(x) / n_pts
        y_mean = sum(y) / n_pts

        ss_xx = sum((xi - x_mean) ** 2 for xi in x)
        ss_xy = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))

        if ss_xx == 0:
            raise ValueError("Shear rate values must not be all identical.")

        n_fit = ss_xy / ss_xx
        ln_K = y_mean - n_fit * x_mean
        K_fit = math.exp(ln_K)

        # Standard error of n and K in log space
        residuals_log = [yi - (ln_K + n_fit * xi) for xi, yi in zip(x, y)]
        s2_log = sum(r ** 2 for r in residuals_log) / max(1, n_pts - 2)
        se_n = math.sqrt(s2_log / ss_xx) if ss_xx > 0 else 0.0

        # Goodness-of-fit in real stress space
        y_real = [pt[1] for pt in points]
        y_pred = [K_fit * (pt[0] ** n_fit) for pt in points]
        mean_y_real = sum(y_real) / n_pts
        ss_tot = sum((yi - mean_y_real) ** 2 for yi in y_real)
        ss_res = sum((yi - yhat) ** 2 for yi, yhat in zip(y_real, y_pred))

        r_squared = max(0.0, 1.0 - (ss_res / ss_tot)) if ss_tot > 0 else 1.0
        rmse = math.sqrt(ss_res / n_pts)

        return {
            "model": "Ostwald-de Waele (Power-Law)",
            "K": round(K_fit, 4),
            "n": round(n_fit, 4),
            "se_n": round(se_n, 4),
            "yield_stress_pa": 0.0,
            "r_squared": round(r_squared, 4),
            "rmse": round(rmse, 3),
            "points_count": n_pts,
        }

    @staticmethod
    def fit_herschel_bulkley(points):
        """
        Fits Herschel-Bulkley model: tau = tau_y + K * gamma_dot^n
        Finds optimal tau_y in [0, 0.95 * min(tau)] that minimizes real-space sum of squared errors.
        """
        if len(points) < 4:
            # Fall back to power law if insufficient points for 3-parameter model
            return RheologyCalibrator.fit_power_law(points)

        min_tau = min(pt[1] for pt in points)
        tau_y_max = max(0.0, min_tau * 0.95)

        best_fit = None
        best_ss_res = 1e12

        # 100-step grid scan across possible yield stress bounds
        steps = 100
        for i in range(steps + 1):
            tau_y_test = (tau_y_max * i) / steps
            shifted_points = [(pt[0], pt[1] - tau_y_test) for pt in points if pt[1] > tau_y_test]
            if len(shifted_points) < 3:
                continue

            try:
                fit = RheologyCalibrator.fit_power_law(shifted_points)
                # Compute real space residuals: tau_i - (tau_y + K * gamma^n)
                y_pred = [tau_y_test + fit["K"] * (pt[0] ** fit["n"]) for pt in points]
                ss_res = sum((pt[1] - yhat) ** 2 for pt, yhat in zip(points, y_pred))

                if ss_res < best_ss_res:
                    best_ss_res = ss_res
                    best_fit = {
                        "model": "Herschel-Bulkley",
                        "yield_stress_pa": round(tau_y_test, 3),
                        "K": fit["K"],
                        "n": fit["n"],
                        "se_n": fit["se_n"],
                        "ss_res": ss_res,
                    }
            except Exception:
                continue

        if not best_fit:
            return RheologyCalibrator.fit_power_law(points)

        # Calculate final goodness-of-fit
        n_pts = len(points)
        y_real = [pt[1] for pt in points]
        mean_y = sum(y_real) / n_pts
        ss_tot = sum((y - mean_y) ** 2 for y in y_real)
        r_squared = max(0.0, 1.0 - (best_fit["ss_res"] / ss_tot)) if ss_tot > 0 else 1.0
        rmse = math.sqrt(best_fit["ss_res"] / n_pts)

        best_fit["r_squared"] = round(r_squared, 4)
        best_fit["rmse"] = round(rmse, 3)
        best_fit["points_count"] = n_pts
        del best_fit["ss_res"]

        return best_fit

    @staticmethod
    def generate_comparison_curve(points, fit_result):
        """
        Generates 50 smooth synthetic points for plotting the fitted model against observed data.
        """
        min_gamma = min(pt[0] for pt in points)
        max_gamma = max(pt[0] for pt in points)
        K = fit_result["K"]
        n = fit_result["n"]
        tau_y = fit_result.get("yield_stress_pa", 0.0)

        # Logarithmic spacing for smooth visualization across decades
        log_min = math.log10(max(0.01, min_gamma))
        log_max = math.log10(max_gamma)
        num_curve_pts = 40
        step = (log_max - log_min) / (num_curve_pts - 1) if num_curve_pts > 1 else 1.0

        curve = []
        for i in range(num_curve_pts):
            gamma = 10.0 ** (log_min + i * step)
            pred_stress = tau_y + K * (gamma ** n)
            pred_visc = pred_stress / gamma
            curve.append({
                "shear_rate": round(gamma, 2),
                "shear_stress_pa": round(pred_stress, 3),
                "viscosity_pa_s": round(pred_visc, 4),
            })

        return curve

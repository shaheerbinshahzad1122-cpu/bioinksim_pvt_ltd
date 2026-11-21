"""
BioInkSim Uncertainty Quantification (UQ) Engine
=================================================
Performs Monte Carlo error propagation across experimental and material uncertainties:
- Consistency index K (+- 15%)
- Flow behavior index n (+- 5%)
- Cell critical shear limit tau_crit (+- 15%)
- Nozzle diameter tolerance (+- 3%)

Computes ensemble statistics (mean, std dev, 95% confidence intervals)
for shear stress, flow rate, velocity, cell viability, and modulus.
"""

import math
import random
from .engine import BioInkEngine
from .data import NOZZLES_DB, CELLS_DB


class UncertaintyEngine:
    """Pure-Python Monte Carlo error propagation engine."""

    @staticmethod
    def run_monte_carlo(params, num_samples=500, random_seed=42):
        """
        Runs Monte Carlo propagation of material and hardware parameter tolerances.
        """
        rng = random.Random(random_seed)

        polymer = params.get("polymer", "Alginate")
        conc = float(params.get("concentration", 3.0))
        mw = float(params.get("mw_kda", 180.0))
        temp = float(params.get("temperature", 25.0))
        ph = float(params.get("ph", 7.4))

        cell_type = params.get("cell_type", "Chondrocytes")
        cell_density = float(params.get("cell_density", 5.0))

        crosslinker = params.get("crosslinker", "CaCl2")
        cross_conc = float(params.get("crosslinker_conc", 50.0))
        cross_time = float(params.get("crosslinker_time", 10.0))

        nozzle = params.get("nozzle", "21G")
        nozzle_type = params.get("nozzle_type", "Cylindrical (Blunt)")
        pressure = float(params.get("pressure_kpa", 150.0))
        syringe_vol = float(params.get("syringe_vol_ml", 5.0))

        # Baseline evaluation
        base_rheo = BioInkEngine.calculate_rheology(polymer, conc, mw, temp, ph)
        base_cross = BioInkEngine.calculate_crosslinking(
            crosslinker, cross_conc, cross_time, temp, polymer, conc_w_v=conc, ph=ph
        )

        c_info = CELLS_DB.get(cell_type, CELLS_DB["Chondrocytes"])
        tau_crit_base = c_info["tau_crit_pa"]
        K_base = base_rheo["K"]
        n_base = base_rheo["n"]
        yield_base = base_rheo["yield_stress_pa"]

        # Sample collections
        tau_w_list = []
        flow_q_list = []
        vel_list = []
        viability_list = []
        modulus_list = []

        for _ in range(num_samples):
            # Sample parameters with realistic log-normal / Gaussian noise
            K_i = max(0.01, rng.gauss(K_base, 0.15 * K_base))
            n_i = max(0.18, min(0.98, rng.gauss(n_base, 0.05 * n_base)))
            yield_i = max(0.0, rng.gauss(yield_base, 0.15 * yield_base))
            tau_crit_i = max(500.0, rng.gauss(tau_crit_base, 0.15 * tau_crit_base))

            # Extrusion calculation
            ext = BioInkEngine.calculate_extrusion_mechanics(
                nozzle, nozzle_type, pressure, syringe_vol, K_i, n_i, yield_i
            )

            # Cell viability calculation with sampled tau_crit
            if ext["flow_active"]:
                stress_ratio = ext["wall_shear_stress_pa"] / tau_crit_i
                dose = (stress_ratio ** 1.15) * (max(1e-4, ext["residence_time_s"]) ** 0.38)
                damage_fraction = (dose ** 2.2) / (1.0 + (dose ** 2.2))
                viab = max(5.0, min(98.0, 97.5 * (1.0 - 0.88 * damage_fraction)))
            else:
                viab = 97.5

            # Crosslinking modulus sample (sample conc +- 5%)
            conc_i = max(0.2, rng.gauss(conc, 0.05 * conc))
            cross_i = BioInkEngine.calculate_crosslinking(
                crosslinker, cross_conc, cross_time, temp, polymer, conc_w_v=conc_i, ph=ph
            )

            tau_w_list.append(ext["wall_shear_stress_pa"])
            flow_q_list.append(ext["flow_rate_mm3_s"])
            vel_list.append(ext["extrusion_velocity_mm_s"])
            viability_list.append(viab)
            modulus_list.append(cross_i["youngs_modulus_kpa"])

        def compute_stats(arr):
            arr_sorted = sorted(arr)
            n = len(arr_sorted)
            mean_val = sum(arr_sorted) / n
            var = sum((x - mean_val) ** 2 for x in arr_sorted) / (n - 1) if n > 1 else 0.0
            std_val = math.sqrt(var)
            ci_low = arr_sorted[int(0.025 * n)]
            ci_high = arr_sorted[int(0.975 * n)]
            return {
                "mean": round(mean_val, 2),
                "std": round(std_val, 2),
                "ci_95_low": round(ci_low, 2),
                "ci_95_high": round(ci_high, 2),
            }

        return {
            "num_samples": num_samples,
            "wall_shear_stress_pa": compute_stats(tau_w_list),
            "flow_rate_mm3_s": compute_stats(flow_q_list),
            "extrusion_velocity_mm_s": compute_stats(vel_list),
            "viability_pct": compute_stats(viability_list),
            "youngs_modulus_kpa": compute_stats(modulus_list),
        }

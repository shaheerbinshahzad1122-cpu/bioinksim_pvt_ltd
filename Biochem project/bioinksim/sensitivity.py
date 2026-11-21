"""
BioInkSim Sensitivity Analysis Engine
=====================================
Performs One-at-a-Time (OAT) parameter sensitivity scans (+- 20%)
to construct interactive Tornado plots showing which physical knobs
drive the greatest change in Composite Score, Viability, and Modulus.
"""

from .engine import BioInkEngine
from .data import NOZZLES_DB


class SensitivityEngine:
    """One-at-a-Time (OAT) local sensitivity scanner."""

    @staticmethod
    def run_sensitivity_scan(base_params):
        """
        Scans parameters +-20% and measures swings in composite score and viability.
        """
        # Baseline run
        def eval_point(p):
            poly = p.get("polymer", "Alginate")
            c = float(p.get("concentration", 3.0))
            mw = float(p.get("mw_kda", 180.0))
            temp = float(p.get("temperature", 25.0))
            ph = float(p.get("ph", 7.4))

            cell_type = p.get("cell_type", "Chondrocytes")
            cell_density = float(p.get("cell_density", 5.0))
            gf = p.get("growth_factor", "TGF-β")

            crosslinker = p.get("crosslinker", "CaCl2")
            cross_conc = float(p.get("crosslinker_conc", 50.0))
            cross_time = float(p.get("crosslinker_time", 10.0))

            nozzle = p.get("nozzle", "21G")
            nozzle_type = p.get("nozzle_type", "Cylindrical (Blunt)")
            pressure = float(p.get("pressure_kpa", 150.0))
            syringe_vol = float(p.get("syringe_vol_ml", 5.0))

            target_tissue = p.get("target_tissue", "Cartilage")
            layer_height = float(p.get("layer_height_um", 250.0))
            strand_spacing = float(p.get("strand_spacing_um", 1000.0))

            rheo = BioInkEngine.calculate_rheology(poly, c, mw, temp, ph)
            cross = BioInkEngine.calculate_crosslinking(
                crosslinker, cross_conc, cross_time, temp, poly, conc_w_v=c, ph=ph
            )
            ext = BioInkEngine.calculate_extrusion_mechanics(
                nozzle, nozzle_type, pressure, syringe_vol, rheo["K"], rheo["n"], rheo["yield_stress_pa"]
            )
            safety = BioInkEngine.calculate_cell_safety(
                cell_type, cell_density, ext["wall_shear_stress_pa"], ext["residence_time_s"], ph, ext["flow_active"]
            )
            scaff = BioInkEngine.calculate_scaffold_metrics(
                target_tissue, nozzle, layer_height, strand_spacing, cell_density,
                rheo["n"], rheo["yield_stress_pa"], cross["crosslinking_degree_pct"]
            )
            deg = BioInkEngine.calculate_degradation_and_swelling(poly, cross["crosslinking_degree_pct"], 37.0, ph)
            suit = BioInkEngine.evaluate_tissue_suitability(
                target_tissue, poly, cell_type, gf, rheo, safety, scaff, cross, deg
            )

            return {
                "score": suit["total_score"],
                "viability": safety["viability_pct"],
                "modulus": cross["youngs_modulus_kpa"],
                "flow_rate": ext["flow_rate_mm3_s"],
            }

        base_res = eval_point(base_params)
        base_score = base_res["score"]

        # Parameters to perturb
        perturbations = [
            {
                "id": "pressure",
                "label": "Extrusion Pressure (kPa)",
                "key": "pressure_kpa",
                "default": float(base_params.get("pressure_kpa", 150.0)),
                "delta_pct": 0.20,
            },
            {
                "id": "concentration",
                "label": "Polymer Concentration (% w/v)",
                "key": "concentration",
                "default": float(base_params.get("concentration", 3.0)),
                "delta_pct": 0.20,
            },
            {
                "id": "crosslinker_conc",
                "label": "Crosslinker Concentration (mM)",
                "key": "crosslinker_conc",
                "default": float(base_params.get("crosslinker_conc", 50.0)),
                "delta_pct": 0.25,
            },
            {
                "id": "crosslinker_time",
                "label": "Crosslinking Duration (min)",
                "key": "crosslinker_time",
                "default": float(base_params.get("crosslinker_time", 10.0)),
                "delta_pct": 0.30,
            },
            {
                "id": "temperature",
                "label": "Temperature (°C)",
                "key": "temperature",
                "default": float(base_params.get("temperature", 25.0)),
                "delta_pct": 0.15,
            },
            {
                "id": "mw",
                "label": "Molecular Weight (kDa)",
                "key": "mw_kda",
                "default": float(base_params.get("mw_kda", 180.0)),
                "delta_pct": 0.20,
            },
        ]

        tornado_items = []
        for item in perturbations:
            val_base = item["default"]
            delta = item["delta_pct"]

            val_low = max(0.1, val_base * (1.0 - delta))
            val_high = val_base * (1.0 + delta)

            p_low = dict(base_params)
            p_low[item["key"]] = val_low
            res_low = eval_point(p_low)

            p_high = dict(base_params)
            p_high[item["key"]] = val_high
            res_high = eval_point(p_high)

            swing_score = abs(res_high["score"] - res_low["score"])
            tornado_items.append({
                "id": item["id"],
                "label": item["label"],
                "base_val": round(val_base, 2),
                "low_val": round(val_low, 2),
                "high_val": round(val_high, 2),
                "score_low": res_low["score"],
                "score_high": res_high["score"],
                "delta_score_low": round(res_low["score"] - base_score, 1),
                "delta_score_high": round(res_high["score"] - base_score, 1),
                "viability_low": res_low["viability"],
                "viability_high": res_high["viability"],
                "modulus_low": res_low["modulus"],
                "modulus_high": res_high["modulus"],
                "swing_score": round(swing_score, 2),
            })

        # Sort descending by swing in score for classic tornado plot presentation
        tornado_items.sort(key=lambda x: x["swing_score"], reverse=True)

        return {
            "baseline_score": base_score,
            "baseline_viability": base_res["viability"],
            "baseline_modulus": base_res["modulus"],
            "tornado_items": tornado_items,
        }

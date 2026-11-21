"""
BioInkSim HTTP REST API & Web Application Server
================================================
Multi-threaded HTTP server providing REST API endpoints and serving frontend assets
(HTML, CSS, JS) with robust error handling, input validation, and security headers.
"""

import http.server
import socketserver
import json
import os
import urllib.parse
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


FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "frontend")


class BioInkRequestHandler(http.server.SimpleHTTPRequestHandler):
    """Custom request handler serving REST endpoints and single-page application assets."""

    def end_headers(self):
        # Security & caching headers
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "SAMEORIGIN")
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        super().end_headers()

    def send_json(self, data, status_code=200):
        body = json.dumps(data, indent=2).encode("utf-8")
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_error_json(self, status_code, message, param=None):
        resp = {
            "status": "error",
            "code": status_code,
            "message": message,
            "param": param,
        }
        self.send_json(resp, status_code=status_code)

    def serve_static_file(self, filename, content_type):
        filepath = os.path.join(FRONTEND_DIR, filename)
        if not os.path.exists(filepath):
            self.send_error(404, f"File {filename} not found.")
            return

        with open(filepath, "rb") as f:
            content = f.read()

        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def _compute_full_simulation(self, req_data):
        """Helper to run a complete forward simulation pipeline from parameters."""
        polymer = req_data.get("polymer", "Alginate")
        conc = float(req_data.get("concentration", 3.0))
        mw = float(req_data.get("mw_kda", 180.0))
        temp = float(req_data.get("temperature", 25.0))
        ph = float(req_data.get("ph", 7.4))

        cell_type = req_data.get("cell_type", "Chondrocytes")
        cell_density = float(req_data.get("cell_density", 5.0))
        gf = req_data.get("growth_factor", "TGF-β")

        crosslinker = req_data.get("crosslinker", "CaCl2")
        cross_conc = float(req_data.get("crosslinker_conc", 50.0))
        cross_time = float(req_data.get("crosslinker_time", 10.0))

        nozzle = req_data.get("nozzle", "21G")
        nozzle_type = req_data.get("nozzle_type", "Cylindrical (Blunt)")
        pressure = float(req_data.get("pressure_kpa", 150.0))
        syringe_vol = float(req_data.get("syringe_vol_ml", 5.0))

        target_tissue = req_data.get("target_tissue", "Cartilage")
        layer_height = float(req_data.get("layer_height_um", 250.0))
        strand_spacing = float(req_data.get("strand_spacing_um", 1000.0))

        rheo = BioInkEngine.calculate_rheology(polymer, conc, mw, temp, ph)
        cross = BioInkEngine.calculate_crosslinking(
            crosslinker, cross_conc, cross_time, temp, polymer, conc_w_v=conc, ph=ph
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
        deg = BioInkEngine.calculate_degradation_and_swelling(
            polymer, cross["crosslinking_degree_pct"], 37.0, ph
        )
        eval_data = BioInkEngine.evaluate_tissue_suitability(
            target_tissue, polymer, cell_type, gf, rheo, safety, scaff, cross, deg
        )

        return {
            "inputs": {
                "polymer": polymer,
                "concentration": conc,
                "mw_kda": mw,
                "temperature": temp,
                "ph": ph,
                "cell_type": cell_type,
                "cell_density": cell_density,
                "growth_factor": gf,
                "crosslinker": crosslinker,
                "crosslinker_conc": cross_conc,
                "crosslinker_time": cross_time,
                "nozzle": nozzle,
                "nozzle_type": nozzle_type,
                "pressure_kpa": pressure,
                "syringe_vol_ml": syringe_vol,
                "target_tissue": target_tissue,
                "layer_height_um": layer_height,
                "strand_spacing_um": strand_spacing,
            },
            "rheology": rheo,
            "crosslinking": cross,
            "extrusion": ext,
            "cell_safety": safety,
            "scaffold": scaff,
            "degradation": deg,
            "tissue_evaluation": eval_data,
        }

    def do_GET(self):
        url_parts = urllib.parse.urlparse(self.path)
        path = url_parts.path

        if path in ["/", "/index.html"]:
            self.serve_static_file("index.html", "text/html; charset=utf-8")
            return

        elif path == "/style.css":
            self.serve_static_file("style.css", "text/css; charset=utf-8")
            return

        elif path == "/app.js":
            self.serve_static_file("app.js", "application/javascript; charset=utf-8")
            return

        elif path == "/api/presets":
            data = {
                "polymers": list(POLYMERS_DB.keys()),
                "cells": list(CELLS_DB.keys()),
                "growth_factors": list(GROWTH_FACTORS_DB.keys()),
                "nozzles": list(NOZZLES_DB.keys()),
                "crosslinkers": list(CROSSLINKERS_DB.keys()),
                "tissues": list(TISSUES_DB.keys()),
                "polymers_detail": POLYMERS_DB,
                "tissues_detail": TISSUES_DB,
                "benchmarks": BENCHMARKS_DB,
            }
            self.send_json(data)
            return

        elif path == "/api/report/pdf":
            # Direct GET request allows instant download of PDF report in browser
            sim_data = self._compute_full_simulation({})
            pdf_bytes = generate_pdf_report(sim_data)
            self.send_response(200)
            self.send_header("Content-Type", "application/pdf")
            self.send_header("Content-Disposition", 'attachment; filename="BioInkSim_Laboratory_Report.pdf"')
            self.send_header("Content-Length", str(len(pdf_bytes)))
            self.end_headers()
            self.wfile.write(pdf_bytes)
            return

        elif path == "/api/benchmarks":
            # Return benchmarks with live simulator parity comparisons
            benchmarks_with_pred = []
            for b in BENCHMARKS_DB:
                inp = b["inputs"]
                try:
                    rheo = BioInkEngine.calculate_rheology(
                        inp["polymer"], inp["concentration"], inp["mw_kda"], inp["temperature"], inp["ph"]
                    )
                    cross = BioInkEngine.calculate_crosslinking(
                        inp["crosslinker"], inp["crosslinker_conc"], inp["crosslinker_time"],
                        inp["temperature"], inp["polymer"], conc_w_v=inp["concentration"], ph=inp["ph"]
                    )
                    ext = BioInkEngine.calculate_extrusion_mechanics(
                        inp["nozzle"], inp["nozzle_type"], inp["pressure_kpa"], 5.0,
                        rheo["K"], rheo["n"], rheo["yield_stress_pa"]
                    )
                    safety = BioInkEngine.calculate_cell_safety(
                        inp["cell_type"], 5.0, ext["wall_shear_stress_pa"], ext["residence_time_s"], inp["ph"], ext["flow_active"]
                    )
                    scaff = BioInkEngine.calculate_scaffold_metrics(
                        "Cartilage", inp["nozzle"], 250.0, 1000.0, 5.0, rheo["n"], rheo["yield_stress_pa"], cross["crosslinking_degree_pct"]
                    )
                    deg = BioInkEngine.calculate_degradation_and_swelling(inp["polymer"], cross["crosslinking_degree_pct"], 37.0, inp["ph"])

                    predicted = {
                        "K_pa_s_n": rheo["K"],
                        "n": rheo["n"],
                        "yield_stress_pa": rheo["yield_stress_pa"],
                        "strand_diameter_mm": round(scaff["strand_diameter_um"] / 1000.0, 2),
                        "printability_pr": scaff["printability_pr"],
                        "youngs_modulus_kpa": cross["youngs_modulus_kpa"],
                        "viability_pct": safety["viability_pct"],
                        "half_life_days": deg["half_life_days"],
                    }
                except Exception as ex:
                    predicted = {"error": str(ex)}

                benchmarks_with_pred.append({
                    "id": b["id"],
                    "title": b["title"],
                    "formulation": b["formulation"],
                    "inputs": b["inputs"],
                    "measured": b["measured"],
                    "predicted": predicted,
                    "notes": b["notes"],
                })

            self.send_json(benchmarks_with_pred)
            return

        else:
            self.send_error(404, "Endpoint Not Found")

    def do_POST(self):
        url_parts = urllib.parse.urlparse(self.path)
        path = url_parts.path

        content_length = int(self.headers.get("Content-Length", 0))
        post_body = self.rfile.read(content_length).decode("utf-8") if content_length > 0 else "{}"

        try:
            req_data = json.loads(post_body) if post_body.strip() else {}
        except Exception as e:
            self.send_error_json(400, f"Malformed JSON request body: {str(e)}")
            return

        if path == "/api/report/pdf":
            try:
                # If full calculated simulation is provided, use it; otherwise compute from inputs
                if "inputs" in req_data and "rheology" in req_data and "crosslinking" in req_data:
                    sim_data = req_data
                else:
                    sim_data = self._compute_full_simulation(req_data)

                pdf_bytes = generate_pdf_report(sim_data)
                polymer = sim_data.get("inputs", {}).get("polymer", "Bioink").replace(" ", "_")
                tissue = sim_data.get("inputs", {}).get("target_tissue", "Construct").replace(" ", "_")
                filename = f"BioInkSim_Laboratory_Report_{polymer}_{tissue}.pdf"

                self.send_response(200)
                self.send_header("Content-Type", "application/pdf")
                self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
                self.send_header("Content-Length", str(len(pdf_bytes)))
                self.end_headers()
                self.wfile.write(pdf_bytes)
                return
            except Exception as e:
                self.send_error_json(400, f"PDF report generation failed: {str(e)}")
                return

        elif path == "/api/calculate":
            try:
                polymer = req_data.get("polymer", "Alginate")
                conc = float(req_data.get("concentration", 3.0))
                mw = float(req_data.get("mw_kda", 180.0))
                temp = float(req_data.get("temperature", 25.0))
                ph = float(req_data.get("ph", 7.4))

                cell_type = req_data.get("cell_type", "Chondrocytes")
                cell_density = float(req_data.get("cell_density", 5.0))
                gf = req_data.get("growth_factor", "TGF-β")

                crosslinker = req_data.get("crosslinker", "CaCl2")
                cross_conc = float(req_data.get("crosslinker_conc", 50.0))
                cross_time = float(req_data.get("crosslinker_time", 10.0))

                nozzle = req_data.get("nozzle", "21G")
                nozzle_type = req_data.get("nozzle_type", "Cylindrical (Blunt)")
                pressure = float(req_data.get("pressure_kpa", 150.0))
                syringe_vol = float(req_data.get("syringe_vol_ml", 5.0))

                target_tissue = req_data.get("target_tissue", "Cartilage")
                layer_height = float(req_data.get("layer_height_um", 250.0))
                strand_spacing = float(req_data.get("strand_spacing_um", 1000.0))

                # Bounds checks
                if conc <= 0 or mw <= 0 or pressure <= 0:
                    raise ValueError("Concentration, molecular weight, and extrusion pressure must be positive.")

                # Calculations
                rheo = BioInkEngine.calculate_rheology(polymer, conc, mw, temp, ph)
                cross = BioInkEngine.calculate_crosslinking(
                    crosslinker, cross_conc, cross_time, temp, polymer, conc_w_v=conc, ph=ph
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
                deg = BioInkEngine.calculate_degradation_and_swelling(
                    polymer, cross["crosslinking_degree_pct"], 37.0, ph
                )
                eval_data = BioInkEngine.evaluate_tissue_suitability(
                    target_tissue, polymer, cell_type, gf, rheo, safety, scaff, cross, deg
                )

                # Diagnostic warnings
                diagnostics = BioInkEngine.generate_diagnostics(
                    polymer, crosslinker, rheo, ext, safety, scaff, cross, target_tissue, ph
                )

                # Viscosity curve data points (0.5 to 1500 s^-1)
                curve_points = []
                for gamma in [0.5, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 400.0, 750.0, 1000.0, 1500.0]:
                    stress = rheo["yield_stress_pa"] + rheo["K"] * (gamma ** rheo["n"])
                    visc = stress / gamma
                    curve_points.append({
                        "shear_rate": gamma,
                        "shear_stress_pa": round(stress, 2),
                        "viscosity_pa_s": round(visc, 3),
                    })

                resp = {
                    "status": "success",
                    "rheology": rheo,
                    "crosslinking": cross,
                    "extrusion": ext,
                    "cell_safety": safety,
                    "scaffold": scaff,
                    "degradation": deg,
                    "tissue_evaluation": eval_data,
                    "diagnostics": diagnostics,
                    "rheology_curve": curve_points,
                }
                self.send_json(resp)
                return

            except Exception as e:
                self.send_error_json(400, f"Simulation calculation failed: {str(e)}")
                return

        elif path == "/api/calibrate":
            try:
                raw_data = req_data.get("raw_data", "")
                data_type = req_data.get("data_type", "stress")

                points = RheologyCalibrator.parse_data(raw_data, data_type=data_type)
                if len(points) < 3:
                    self.send_error_json(400, "At least 3 valid (shear_rate, shear_stress/viscosity) points are required.")
                    return

                power_law_fit = RheologyCalibrator.fit_power_law(points)
                hb_fit = RheologyCalibrator.fit_herschel_bulkley(points)
                chosen_fit = hb_fit if (hb_fit["r_squared"] >= power_law_fit["r_squared"] - 0.01 and hb_fit.get("yield_stress_pa", 0) > 0.05) else power_law_fit

                comparison_curve = RheologyCalibrator.generate_comparison_curve(points, chosen_fit)
                exp_points = [{"shear_rate": pt[0], "shear_stress_pa": round(pt[1], 3), "viscosity_pa_s": round(pt[1]/pt[0], 4)} for pt in points]

                resp = {
                    "status": "success",
                    "best_fit": chosen_fit,
                    "power_law_fit": power_law_fit,
                    "herschel_bulkley_fit": hb_fit,
                    "experimental_points": exp_points,
                    "fitted_curve": comparison_curve,
                }
                self.send_json(resp)
                return

            except Exception as e:
                self.send_error_json(400, f"Calibration failed: {str(e)}")
                return

        elif path == "/api/uq":
            try:
                samples = int(req_data.get("samples", 400))
                uq_result = UncertaintyEngine.run_monte_carlo(req_data, num_samples=min(1000, max(100, samples)))
                self.send_json({"status": "success", "uq": uq_result})
                return
            except Exception as e:
                self.send_error_json(400, f"Uncertainty quantification failed: {str(e)}")
                return

        elif path == "/api/sensitivity":
            try:
                sens_result = SensitivityEngine.run_sensitivity_scan(req_data)
                self.send_json({"status": "success", "sensitivity": sens_result})
                return
            except Exception as e:
                self.send_error_json(400, f"Sensitivity analysis failed: {str(e)}")
                return

        elif path == "/api/inverse_design":
            try:
                target_tissue = req_data.get("target_tissue", "Cartilage")
                desired_porosity = float(req_data.get("desired_porosity", 70.0))
                desired_modulus = float(req_data.get("desired_modulus", 180.0))
                layer_height = float(req_data.get("layer_height_um", 250.0))

                optimal = BioInkEngine.run_inverse_design(
                    target_tissue, desired_porosity, desired_modulus, layer_height
                )
                if optimal:
                    self.send_json({"status": "success", "optimal": optimal})
                else:
                    self.send_error_json(404, "No viable formulation found for given targets.")
                return

            except Exception as e:
                self.send_error_json(400, f"Inverse design failed: {str(e)}")
                return

        else:
            self.send_error(404, "Endpoint Not Found")

    def log_message(self, format, *args):
        # Suppress noisy standard HTTP access logs
        pass


class ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True

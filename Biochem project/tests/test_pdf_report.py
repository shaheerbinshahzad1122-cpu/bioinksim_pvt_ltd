"""
Tests for PDF Laboratory Report Generation
==========================================
"""

from bioinksim.engine import BioInkEngine
from bioinksim.report_pdf import generate_pdf_report


def test_pdf_report_compilation():
    """Verify that generate_pdf_report compiles a valid multi-page PDF document."""
    params = {
        "polymer": "Alginate", "concentration": 3.0, "mw_kda": 180, "temperature": 25.0, "ph": 7.4,
        "cell_type": "Chondrocytes", "cell_density": 5.0, "growth_factor": "TGF-β",
        "crosslinker": "CaCl2", "crosslinker_conc": 50.0, "crosslinker_time": 10.0,
        "nozzle": "21G", "nozzle_type": "Cylindrical (Blunt)", "pressure_kpa": 150.0, "syringe_vol_ml": 5.0,
        "target_tissue": "Cartilage", "layer_height_um": 250.0, "strand_spacing_um": 1000.0
    }
    rheo = BioInkEngine.calculate_rheology(params["polymer"], params["concentration"], params["mw_kda"], params["temperature"], params["ph"])
    cross = BioInkEngine.calculate_crosslinking(params["crosslinker"], params["crosslinker_conc"], params["crosslinker_time"], params["temperature"], params["polymer"], conc_w_v=params["concentration"], ph=params["ph"])
    ext = BioInkEngine.calculate_extrusion_mechanics(params["nozzle"], params["nozzle_type"], params["pressure_kpa"], params["syringe_vol_ml"], rheo["K"], rheo["n"], rheo["yield_stress_pa"])
    safety = BioInkEngine.calculate_cell_safety(params["cell_type"], params["cell_density"], ext["wall_shear_stress_pa"], ext["residence_time_s"], params["ph"], ext["flow_active"])
    scaff = BioInkEngine.calculate_scaffold_metrics(params["target_tissue"], params["nozzle"], params["layer_height_um"], params["strand_spacing_um"], params["cell_density"], rheo["n"], rheo["yield_stress_pa"], cross["crosslinking_degree_pct"])
    deg = BioInkEngine.calculate_degradation_and_swelling(params["polymer"], cross["crosslinking_degree_pct"], 37.0, params["ph"])
    suit = BioInkEngine.evaluate_tissue_suitability(params["target_tissue"], params["polymer"], params["cell_type"], params["growth_factor"], rheo, safety, scaff, cross, deg)

    sim_data = {
        "inputs": params,
        "rheology": rheo,
        "crosslinking": cross,
        "extrusion": ext,
        "cell_safety": safety,
        "scaffold": scaff,
        "degradation": deg,
        "tissue_evaluation": suit,
    }

    pdf_bytes = generate_pdf_report(sim_data)

    # Standard PDF assertions
    assert isinstance(pdf_bytes, bytes), "PDF output must be raw bytes"
    assert pdf_bytes.startswith(b"%PDF-1.4"), "Header must be %PDF-1.4"
    assert b"%%EOF" in pdf_bytes, "PDF must end with EOF marker"
    assert b"/Type /Pages" in pdf_bytes
    assert b"/Count 2" in pdf_bytes, "Must be a 2-page document"
    assert len(pdf_bytes) > 5000, "PDF should contain substantial content"

"""
BioInkSim Knowledge Base & Literature Database
==============================================
Empirical literature data, biopolymer parameters, cell mechanical limits,
growth factors, nozzles, crosslinkers, tissue benchmarks, and validation datasets.

Primary literature citations:
- Bagley, E. B. (1957). End corrections in the capillary flow of polyethylene. J. Appl. Phys., 28(5), 624-627.
- Herschel, W. H., & Bulkley, R. (1926). Konsistenzmessungen von Gummi-Benzollösungen. Kolloid-Z., 39(4), 291-300.
- Leverett, M. C. et al. (1972). Red blood cell damage by shear stress. Biophysical Journal, 12(3), 257-273.
- Nevaril, C. G. et al. (1968). Erythrocyte damage and destruction induced by shearing stress. J. Lab. Clin. Med., 71(5), 784-790.
- Zhang, Z. et al. (2017). Evaluation of cell viability during extrusion bioprinting. Biofabrication, 9(2), 025011.
- Ning, L. et al. (2018). Influence of extrusion process on cell viability and proliferation. Biofabrication, 10(3), 035010.
- Ouyang, L. et al. (2016). Effect of bioink printability on print fidelity and cell viability. Biofabrication, 8(3), 035020.
- Gao, T. et al. (2018). Optimization of GelMA bioink for 3D bioprinting of vascularized tissue. Biomaterials, 178, 148-162.
- Müller, M. et al. (2015). Nanostructured bioinks for 3D bioprinting. Ann. Biomed. Eng., 43(3), 774-786.
"""

POLYMERS_DB = {
    "Alginate": {
        "name": "Sodium Alginate",
        "type": "Natural Polysaccharide (Brown Algae)",
        "mw_default_kda": 180,
        "conc_min": 0.5, "conc_max": 8.0, "conc_default": 3.0,
        "density_g_cm3": 1.02,
        "base_K": 4.2,             # Consistency index at 3% w/v, 25C (Pa·s^n)
        "conc_exp_K": 2.7,         # K scales as c^2.7
        "base_n": 0.52,            # Flow behavior index
        "yield_stress_base": 1.8,  # Pa at 3% w/v
        "chi_c": 0.12,             # Percolation threshold crosslinking fraction (gel point)
        "degradation_rate_base": 0.012,  # day^-1 (~58 day half-life in physiological saline)
        "k_swell": 0.42,           # day^-1 swelling rate constant
        "swelling_equilibrium": 650,     # % water uptake
        "default_crosslinker": "CaCl2",
        "crosslinking_type": "Ionic",
        "biocompatibility": 95,
        "cell_adhesion": "Low (requires RGD functionalization for integrin binding)",
        "description": "Biocompatible linear copolymer of beta-D-mannuronate (M) and alpha-L-guluronate (G) residues. Undergoes instantaneous ionic gelation in presence of divalent cations (Ca2+) via egg-box chelation."
    },
    "Gelatin": {
        "name": "Gelatin (Type A / Porcine)",
        "type": "Natural Protein (Denatured Collagen)",
        "mw_default_kda": 100,
        "conc_min": 2.0, "conc_max": 15.0, "conc_default": 5.0,
        "density_g_cm3": 1.04,
        "base_K": 18.5,            # At 22C (physical gel state)
        "conc_exp_K": 3.0,
        "base_n": 0.42,
        "yield_stress_base": 32.0,
        "chi_c": 0.15,
        "degradation_rate_base": 0.065,  # Rapid enzymatic & thermal dissolution at 37C
        "k_swell": 0.85,
        "swelling_equilibrium": 450,
        "default_crosslinker": "Thermal / mTG",
        "crosslinking_type": "Thermal / Enzymatic",
        "biocompatibility": 98,
        "cell_adhesion": "Excellent (abundant native RGD motifs for cell attachment)",
        "description": "Thermo-reversible physical hydrogel. Forms triple-helix junction zones below ~28°C. Rapidly liquefies at physiological body temperature (37°C) unless chemically or enzymatically stabilized."
    },
    "GelMA": {
        "name": "Gelatin Methacryloyl (GelMA)",
        "type": "Semi-Synthetic Functionalized Protein",
        "mw_default_kda": 110,
        "conc_min": 3.0, "conc_max": 20.0, "conc_default": 7.0,
        "density_g_cm3": 1.05,
        "base_K": 12.0,
        "conc_exp_K": 2.5,
        "base_n": 0.48,
        "yield_stress_base": 22.0,
        "chi_c": 0.14,
        "degradation_rate_base": 0.015,  # Tunable by degree of methacrylation (DoF)
        "k_swell": 0.35,
        "swelling_equilibrium": 500,
        "default_crosslinker": "LAP + 405nm",
        "crosslinking_type": "Photo-crosslinking",
        "biocompatibility": 96,
        "cell_adhesion": "High (preserves RGD motifs and MMP degradation cleavage sites)",
        "description": "Methacrylamide-modified gelatin. Combines thermo-reversible shear-thinning extrusion behavior at room temperature with stable covalent radical polymerization upon UV/blue light exposure in the presence of photoinitiators (LAP/Irgacure)."
    },
    "Collagen": {
        "name": "Collagen Type I (Rat Tail / Bovine)",
        "type": "Natural Extracellular Matrix Protein",
        "mw_default_kda": 300,
        "conc_min": 0.2, "conc_max": 4.0, "conc_default": 1.2,
        "density_g_cm3": 1.01,
        "base_K": 14.5,
        "conc_exp_K": 2.2,
        "base_n": 0.36,            # Highly shear-thinning fibrillar network
        "yield_stress_base": 16.0,
        "chi_c": 0.10,
        "degradation_rate_base": 0.024,
        "k_swell": 0.50,
        "swelling_equilibrium": 700,
        "default_crosslinker": "Thermal Gelation",
        "crosslinking_type": "Thermal Self-Assembly",
        "biocompatibility": 100,
        "cell_adhesion": "Superior (gold standard native ECM protein for mammalian cells)",
        "description": "Primary structural protein of mammalian connective tissue. Exhibits fibrillogenesis at physiological pH and 37°C. Forms fibrous network with high shear-thinning but low mechanical stiffness unless chemically crosslinked."
    },
    "Chitosan": {
        "name": "Chitosan (Deacetylated Chitin)",
        "type": "Cationic Polysaccharide",
        "mw_default_kda": 150,
        "conc_min": 1.0, "conc_max": 5.0, "conc_default": 2.5,
        "density_g_cm3": 1.03,
        "base_K": 6.8,
        "conc_exp_K": 2.4,
        "base_n": 0.56,
        "yield_stress_base": 5.5,
        "chi_c": 0.15,
        "degradation_rate_base": 0.018,
        "k_swell": 0.28,
        "swelling_equilibrium": 550,
        "default_crosslinker": "Genipin",
        "crosslinking_type": "Ionic / Covalent",
        "biocompatibility": 92,
        "cell_adhesion": "Moderate (cationic charge promotes cell membrane adhesion)",
        "description": "Poly-cationic biopolymer derived from crustacean chitin. Soluble in mild acidic aqueous solutions (pH < 6.2); displays intrinsic antimicrobial, hemostatic, and wound-healing acceleration properties."
    },
    "Hyaluronic Acid": {
        "name": "Hyaluronic Acid (HA / Sodium Hyaluronate)",
        "type": "Natural Glycosaminoglycan (GAG)",
        "mw_default_kda": 800,
        "conc_min": 0.5, "conc_max": 4.0, "conc_default": 1.5,
        "density_g_cm3": 1.02,
        "base_K": 9.5,
        "conc_exp_K": 2.6,
        "base_n": 0.34,            # Strongly pseudoplastic
        "yield_stress_base": 3.2,
        "chi_c": 0.12,
        "degradation_rate_base": 0.035,  # Rapid enzymatic cleavage by hyaluronidase
        "k_swell": 0.65,
        "swelling_equilibrium": 950,
        "default_crosslinker": "LAP + 405nm",
        "crosslinking_type": "Covalent / Photo",
        "biocompatibility": 97,
        "cell_adhesion": "High CD44 receptor-mediated cellular interaction",
        "description": "High-molecular-weight non-sulfated GAG abundant in synovial fluid, vitreous humor, and articular cartilage ECM. Endows bioinks with extraordinary hydration, viscoelastic lubricity, and shear-thinning."
    },
    "Alginate-Gelatin Blend": {
        "name": "Alginate-Gelatin Composite Bioink",
        "type": "Polysaccharide-Protein Interpenetrating Network (IPN)",
        "mw_default_kda": 150,
        "conc_min": 2.0, "conc_max": 10.0, "conc_default": 4.0,
        "density_g_cm3": 1.03,
        "base_K": 15.0,
        "conc_exp_K": 2.6,
        "base_n": 0.45,
        "yield_stress_base": 24.0,
        "chi_c": 0.13,
        "degradation_rate_base": 0.020,
        "k_swell": 0.55,
        "swelling_equilibrium": 520,
        "default_crosslinker": "CaCl2",
        "crosslinking_type": "Dual (Ionic + Thermal)",
        "biocompatibility": 97,
        "cell_adhesion": "High (gelatin provides RGD motifs, alginate provides structural scaffold)",
        "description": "Gold standard composite bioink formulation. Gelatin provides room-temperature physical gelation, shear-thinning fidelity, and RGD cell adhesion sites, while Alginate allows rapid post-printing ionic locking with Ca2+."
    }
}

CELLS_DB = {
    "Chondrocytes": {
        "name": "Articular Chondrocytes",
        "diameter_um": 14.0,
        "tau_crit_pa": 6200,       # Critical shear stress for 50% membrane rupture (Blaeser 2016)
        "k_apop": 0.06,            # Post-print apoptotic rate constant (day^-1)
        "robustness": "Moderate",
        "preferred_tissues": ["Cartilage"],
        "recommended_gf": "TGF-β",
        "default_density": 5.0,    # x 10^6 cells/mL
        "citation": "Blaeser et al. (Biomaterials 2016), Ning et al. (Biofabrication 2018)",
        "description": "Primary specialized cells of healthy cartilage tissue. Sensitive to hydrodynamic pressure, synthesize Type II collagen and aggrecan matrix."
    },
    "Osteoblasts": {
        "name": "Primary Osteoblasts / MC3T3-E1",
        "diameter_um": 18.0,
        "tau_crit_pa": 7200,
        "k_apop": 0.05,
        "robustness": "High",
        "preferred_tissues": ["Bone"],
        "recommended_gf": "BMP-2",
        "default_density": 4.0,
        "citation": "Zhang et al. (Biofabrication 2017)",
        "description": "Bone-forming cells derived from osteoprogenitor mesenchymal lineages. Express alkaline phosphatase (ALP) and synthesize mineralized osteoid matrix."
    },
    "Fibroblasts": {
        "name": "Dermal Fibroblasts (HDFs)",
        "diameter_um": 15.0,
        "tau_crit_pa": 8000,
        "k_apop": 0.04,
        "robustness": "Very High",
        "preferred_tissues": ["Skin", "Muscle"],
        "recommended_gf": "FGF",
        "default_density": 3.0,
        "citation": "Nevaril et al. (1968), Leverett et al. (1972)",
        "description": "Robust connective tissue cells capable of surviving elevated hydrodynamic shear stress; actively deposit Collagen Type I and Type III."
    },
    "Stem Cells": {
        "name": "Mesenchymal Stem Cells (hMSCs)",
        "diameter_um": 20.0,
        "tau_crit_pa": 3400,       # Delicate membrane, large diameter
        "k_apop": 0.09,
        "robustness": "Sensitive",
        "preferred_tissues": ["Cartilage", "Bone", "Muscle", "Blood vessel"],
        "recommended_gf": "TGF-β",
        "default_density": 5.0,
        "citation": "Blaeser et al. (Biomaterials 2016)",
        "description": "Multipotent progenitor cells capable of tri-lineage differentiation. Have large cell diameters and delicate membranes; require gentle extrusion parameters."
    },
    "Endothelial Cells": {
        "name": "Human Umbilical Vein Endothelial Cells (HUVECs)",
        "diameter_um": 13.0,
        "tau_crit_pa": 4400,
        "k_apop": 0.08,
        "robustness": "Moderate-Sensitive",
        "preferred_tissues": ["Blood vessel"],
        "recommended_gf": "VEGF",
        "default_density": 6.0,
        "citation": "Gao et al. (Biomaterials 2018)",
        "description": "Form endothelial lumen monolayers lining vascular networks. Highly responsive to fluid laminar shear stress signaling; essential for angiogenesis."
    }
}

GROWTH_FACTORS_DB = {
    "TGF-β": {
        "name": "Transforming Growth Factor-Beta 1/3 (TGF-β)",
        "target_tissue": "Cartilage",
        "standard_conc_ng_ml": 10.0,
        "half_life_hours": 24.0,
        "biochemical_role": "Upregulates Sox9 transcription, induces chondrogenesis, and stimulates chondrocytes to synthesize Glycosaminoglycans (GAGs) and Collagen Type II."
    },
    "BMP-2": {
        "name": "Bone Morphogenetic Protein-2 (BMP-2)",
        "target_tissue": "Bone",
        "standard_conc_ng_ml": 100.0,
        "half_life_hours": 18.0,
        "biochemical_role": "Potent osteoinductive cytokine; stimulates Runx2 expression, alkaline phosphatase activity, and osteocalcin deposition for mineralization."
    },
    "VEGF": {
        "name": "Vascular Endothelial Growth Factor (VEGF-A)",
        "target_tissue": "Blood vessel",
        "standard_conc_ng_ml": 25.0,
        "half_life_hours": 12.0,
        "biochemical_role": "Key mediator of vasculogenesis and capillary sprout elongation. Enhances endothelial cell chemotaxis and lumen maturation."
    },
    "FGF": {
        "name": "Basic Fibroblast Growth Factor (bFGF / FGF-2)",
        "target_tissue": "Skin",
        "standard_conc_ng_ml": 20.0,
        "half_life_hours": 8.0,
        "biochemical_role": "Promotes robust fibroblast migration, granulation tissue formation, and broad dermal cellular proliferation."
    },
    "None / Control": {
        "name": "Basal Medium (No Exogenous GF)",
        "target_tissue": "General",
        "standard_conc_ng_ml": 0.0,
        "half_life_hours": 0.0,
        "biochemical_role": "Standard cellular baseline control without targeted lineage differentiation signaling."
    }
}

NOZZLES_DB = {
    "18G": {"inner_diam_mm": 0.838, "outer_diam_mm": 1.270, "color": "Pink"},
    "20G": {"inner_diam_mm": 0.603, "outer_diam_mm": 0.908, "color": "Yellow"},
    "21G": {"inner_diam_mm": 0.514, "outer_diam_mm": 0.819, "color": "Dark Green"},
    "22G": {"inner_diam_mm": 0.413, "outer_diam_mm": 0.718, "color": "Black"},
    "23G": {"inner_diam_mm": 0.337, "outer_diam_mm": 0.641, "color": "Deep Blue"},
    "25G": {"inner_diam_mm": 0.260, "outer_diam_mm": 0.515, "color": "Royal Blue"},
    "27G": {"inner_diam_mm": 0.210, "outer_diam_mm": 0.413, "color": "Medium Grey"},
    "30G": {"inner_diam_mm": 0.159, "outer_diam_mm": 0.311, "color": "Light Yellow"}
}

CROSSLINKERS_DB = {
    "CaCl2": {
        "name": "Calcium Chloride (CaCl2)",
        "type": "Ionic Divalent Salt",
        "applicable_polymers": ["Alginate", "Alginate-Gelatin Blend"],
        "default_conc_mm": 50.0,
        "conc_min_mm": 10.0, "conc_max_mm": 200.0,
        "default_time_min": 10.0,
        "mechanism": "Binds electrostatically to G-blocks of alginate chains, creating cooperative egg-box coordination junction complexes.",
        "max_modulus_kpa": 120.0
    },
    "LAP + 405nm": {
        "name": "Lithium Phenyl-2,4,6-trimethylbenzoylphosphinate (LAP)",
        "type": "Water-Soluble Photoinitiator (Visible Blue 405nm / UV 365nm)",
        "applicable_polymers": ["GelMA", "Hyaluronic Acid"],
        "default_conc_mm": 1.7, # ~0.05% w/v (~1.7 mM)
        "conc_min_mm": 0.3, "conc_max_mm": 5.0,
        "default_time_min": 1.0, # 60 seconds
        "mechanism": "Photolytic alpha-cleavage produces benzoyl and phosphonyl radicals initiating chain-growth radical polymerization of methacryloyl double bonds.",
        "max_modulus_kpa": 180.0
    },
    "Thermal Gelation": {
        "name": "Temperature-Induced Self-Assembly",
        "type": "Physical Entanglement & Helix Formation",
        "applicable_polymers": ["Gelatin", "Collagen", "Alginate-Gelatin Blend"],
        "default_conc_mm": 0.0,
        "conc_min_mm": 0.0, "conc_max_mm": 0.0,
        "default_time_min": 15.0,
        "mechanism": "Hydrogen bonding and hydrophobic collapse induce triple-helical coil-to-helix transitions below gelation transition temperature.",
        "max_modulus_kpa": 45.0
    },
    "mTG": {
        "name": "Microbial Transglutaminase (mTG)",
        "type": "Enzymatic Crosslinker",
        "applicable_polymers": ["Gelatin", "Collagen", "Alginate-Gelatin Blend"],
        "default_conc_mm": 20.0, # U/g
        "conc_min_mm": 5.0, "conc_max_mm": 80.0,
        "default_time_min": 25.0,
        "mechanism": "Catalyzes acyl-transfer reaction between gamma-carboxamide group of glutamine residues and epsilon-amino group of lysine residues, forming isopeptide bonds.",
        "max_modulus_kpa": 140.0
    },
    "Genipin": {
        "name": "Genipin (Gardenia jasminoides extract)",
        "type": "Naturally Derived Covalent Bifunctional Crosslinker",
        "applicable_polymers": ["Chitosan", "Gelatin", "Collagen"],
        "default_conc_mm": 5.0,
        "conc_min_mm": 0.5, "conc_max_mm": 20.0,
        "default_time_min": 45.0,
        "mechanism": "Nucleophilic attack by primary amines onto genipin ring followed by dimerization, yielding stable bluish covalent crosslinks with low cytotoxicity.",
        "max_modulus_kpa": 95.0
    }
}

TISSUES_DB = {
    "Cartilage": {
        "name": "Articular Cartilage",
        "target_modulus_kpa_min": 80.0,
        "target_modulus_kpa_max": 400.0,
        "target_modulus_ideal": 200.0,
        "target_porosity_min": 60.0,
        "target_porosity_max": 80.0,
        "target_pore_size_um": 180.0,
        "preferred_cells": ["Chondrocytes", "Stem Cells"],
        "preferred_gf": ["TGF-β"],
        "preferred_polymers": ["Alginate", "GelMA", "Hyaluronic Acid", "Alginate-Gelatin Blend"],
        "degradation_half_life_ideal_days": 60.0,
        "hydration_target": "High (70-85% water content)",
        "clinical_significance": "Avascular tissue with limited self-repair. Requires stable compressive viscoelasticity, chondrocyte survival, and gradual matrix replacement by neo-synthesized aggrecan."
    },
    "Bone": {
        "name": "Trabecular / Cortical Bone",
        "target_modulus_kpa_min": 500.0,
        "target_modulus_kpa_max": 2500.0,
        "target_modulus_ideal": 1200.0,
        "target_porosity_min": 50.0,
        "target_porosity_max": 75.0,
        "target_pore_size_um": 320.0,
        "preferred_cells": ["Osteoblasts", "Stem Cells"],
        "preferred_gf": ["BMP-2"],
        "preferred_polymers": ["GelMA", "Alginate", "Alginate-Gelatin Blend", "Chitosan"],
        "degradation_half_life_ideal_days": 120.0,
        "hydration_target": "Moderate with high mineral nucleation potential",
        "clinical_significance": "Requires high structural stiffness, osteoinductive signaling (BMP-2), and large interconnected pores (300-400 um) for cell infiltration and osteon vascularization."
    },
    "Skin": {
        "name": "Dermis / Dermo-Epidermal Construct",
        "target_modulus_kpa_min": 10.0,
        "target_modulus_kpa_max": 65.0,
        "target_modulus_ideal": 30.0,
        "target_porosity_min": 70.0,
        "target_porosity_max": 88.0,
        "target_pore_size_um": 120.0,
        "preferred_cells": ["Fibroblasts", "Stem Cells"],
        "preferred_gf": ["FGF", "VEGF"],
        "preferred_polymers": ["Collagen", "Gelatin", "GelMA", "Chitosan", "Alginate-Gelatin Blend"],
        "degradation_half_life_ideal_days": 28.0,
        "hydration_target": "High hydration with flexible pliability",
        "clinical_significance": "Demands compliant, flexible hydrogel matching skin elasticity, fast cellular migration, and rapid angiogenesis to prevent ischemic wound necrosis."
    },
    "Muscle": {
        "name": "Skeletal / Cardiac Muscle Tissue",
        "target_modulus_kpa_min": 15.0,
        "target_modulus_kpa_max": 90.0,
        "target_modulus_ideal": 40.0,
        "target_porosity_min": 60.0,
        "target_porosity_max": 80.0,
        "target_pore_size_um": 150.0,
        "preferred_cells": ["Fibroblasts", "Stem Cells"],
        "preferred_gf": ["FGF"],
        "preferred_polymers": ["GelMA", "Collagen", "Alginate-Gelatin Blend"],
        "degradation_half_life_ideal_days": 45.0,
        "hydration_target": "Elastic and highly compliant matrix",
        "clinical_significance": "Requires aligned anisotropic strand architecture, fatigue resistance under cyclic contractile strain, and cell elongation guidance cues."
    },
    "Blood vessel": {
        "name": "Vascular Graft / Endothelial Channel",
        "target_modulus_kpa_min": 25.0,
        "target_modulus_kpa_max": 120.0,
        "target_modulus_ideal": 60.0,
        "target_porosity_min": 65.0,
        "target_porosity_max": 85.0,
        "target_pore_size_um": 100.0,
        "preferred_cells": ["Endothelial Cells", "Stem Cells"],
        "preferred_gf": ["VEGF"],
        "preferred_polymers": ["GelMA", "Alginate-Gelatin Blend", "Collagen"],
        "degradation_half_life_ideal_days": 40.0,
        "hydration_target": "Smooth, non-thrombogenic lumen surface",
        "clinical_significance": "Demands tubular patency, burst pressure resistance, non-thrombogenic inner wall, and sustained VEGF signaling for intact confluent endothelium formation."
    }
}

BENCHMARKS_DB = [
    {
        "id": "ouyang_2016_alginate_gelatin",
        "title": "Ouyang et al. (Biofabrication 2016)",
        "formulation": "3% w/v Alginate + 4% w/v Gelatin @ 22°C, 21G Needle, 120 kPa",
        "inputs": {
            "polymer": "Alginate-Gelatin Blend",
            "concentration": 4.0,
            "mw_kda": 150.0,
            "temperature": 22.0,
            "ph": 7.4,
            "nozzle": "21G",
            "nozzle_type": "Cylindrical (Blunt)",
            "pressure_kpa": 120.0,
            "crosslinker": "CaCl2",
            "crosslinker_conc": 50.0,
            "crosslinker_time": 10.0,
            "cell_type": "Chondrocytes"
        },
        "measured": {
            "K_pa_s_n": 15.0,
            "n": 0.45,
            "strand_diameter_mm": 0.58,
            "printability_pr": 0.98,
            "viability_pct": 89.0
        },
        "notes": "Classic study establishing the Ouyang square-pore printability index (Pr). Demonstrates high cell survival with balanced gelatin/alginate dual crosslinking."
    },
    {
        "id": "gao_2018_gelma",
        "title": "Gao et al. (Biomaterials 2018)",
        "formulation": "7% w/v GelMA @ 23°C, 22G Needle, 100 kPa, LAP 405nm (60s)",
        "inputs": {
            "polymer": "GelMA",
            "concentration": 7.0,
            "mw_kda": 110.0,
            "temperature": 23.0,
            "ph": 7.4,
            "nozzle": "22G",
            "nozzle_type": "Conical (Tapered)",
            "pressure_kpa": 100.0,
            "crosslinker": "LAP + 405nm",
            "crosslinker_conc": 1.7,
            "crosslinker_time": 1.0,
            "cell_type": "Endothelial Cells"
        },
        "measured": {
            "K_pa_s_n": 12.0,
            "n": 0.48,
            "yield_stress_pa": 22.0,
            "youngs_modulus_kpa": 32.0,
            "viability_pct": 86.5
        },
        "notes": "Vascularized hydrogel construct optimization. High shape fidelity combined with rapid photo-crosslinking preserving endothelial cell viability."
    },
    {
        "id": "muller_2015_alginate",
        "title": "Müller et al. (Ann. Biomed. Eng. 2015)",
        "formulation": "2.5% w/v Alginate @ 25°C, 25G Needle, 160 kPa, 50 mM CaCl2 (10 min)",
        "inputs": {
            "polymer": "Alginate",
            "concentration": 2.5,
            "mw_kda": 180.0,
            "temperature": 25.0,
            "ph": 7.4,
            "nozzle": "25G",
            "nozzle_type": "Cylindrical (Blunt)",
            "pressure_kpa": 160.0,
            "crosslinker": "CaCl2",
            "crosslinker_conc": 50.0,
            "crosslinker_time": 10.0,
            "cell_type": "Chondrocytes"
        },
        "measured": {
            "K_pa_s_n": 3.8,
            "n": 0.53,
            "youngs_modulus_kpa": 48.0,
            "half_life_days": 65.0,
            "viability_pct": 82.0
        },
        "notes": "Chondrocyte encapsulation in nanostructured alginate for cartilage repair with long-term mechanical stability."
    }
]

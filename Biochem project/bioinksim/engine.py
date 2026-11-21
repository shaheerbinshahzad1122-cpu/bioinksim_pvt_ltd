"""
BioInkSim Calculation Engine
============================
Biophysically rigorous simulator modeling non-Newtonian rheology, crosslinking kinetics,
capillary extrusion hydrodynamics, cell mechanical viability, 3D scaffold architecture,
degradation kinetics, tissue suitability scoring, and inverse design.
"""

import math
from .data import (
    POLYMERS_DB,
    CELLS_DB,
    GROWTH_FACTORS_DB,
    NOZZLES_DB,
    CROSSLINKERS_DB,
    TISSUES_DB,
)


class BioInkEngine:
    """Core mathematical and biochemical calculation engine for BioInkSim."""

    @staticmethod
    def calculate_rheology(polymer_name, conc_w_v, mw_kda, temp_c, ph):
        """
        Calculates non-Newtonian Power-Law & Herschel-Bulkley parameters.
        Accounts for polymer concentration scaling, molecular weight, Arrhenius temperature shift,
        and polyelectrolyte pH charge shifts.
        """
        p_info = POLYMERS_DB.get(polymer_name, POLYMERS_DB["Alginate"])
        base_K = p_info["base_K"]
        conc_exp = p_info["conc_exp_K"]
        base_n = p_info["base_n"]
        base_mw = p_info["mw_default_kda"]
        yield_base = p_info["yield_stress_base"]

        # 1. Concentration scaling: K ~ c^exp
        c_ratio = max(0.1, conc_w_v / p_info["conc_default"])
        K_c = base_K * (c_ratio ** conc_exp)

        # 2. Molecular weight scaling (Mark-Houwink / Reptation regime):
        # Viscosity ~ (Mw / Mw_0)^1.4 in semi-dilute entangled hydrogel solutions
        mw_ratio = max(0.1, mw_kda / base_mw)
        K_mw = K_c * (mw_ratio ** 1.4)

        # 3. Temperature modulation:
        if polymer_name in ["Gelatin", "GelMA"]:
            # Thermo-reversible triple helix transition around 27.5 °C
            t_gel = 27.5
            if temp_c <= t_gel:
                # Solid-like gelled network
                K_temp = K_mw * math.exp(0.18 * (t_gel - temp_c))
                flow_n = max(0.30, base_n - 0.015 * (t_gel - temp_c))
                yield_stress = yield_base * (c_ratio ** 2.2) * (1.0 + 0.35 * (t_gel - temp_c))
            else:
                # Disordered coil sol state
                delta_t = temp_c - t_gel
                K_temp = K_mw * math.exp(-0.12 * delta_t)
                flow_n = min(0.95, base_n + 0.03 * delta_t)
                yield_stress = max(0.1, yield_base * math.exp(-0.3 * delta_t))
        else:
            # Standard Arrhenius relationship: K(T) = K_0 * exp[ (Ea/R) * (1/T - 1/T_ref) ]
            T_k = temp_c + 273.15
            T_ref_k = 25.0 + 273.15  # 298.15 K
            # Average activation energy for polysaccharide flow Ea/R ~ 2400 K
            K_temp = K_mw * math.exp(2400.0 * (1.0 / T_k - 1.0 / T_ref_k))
            flow_n = max(0.20, min(0.98, base_n - 0.03 * (conc_w_v - p_info["conc_default"])))
            yield_stress = yield_base * (c_ratio ** 1.8)

        # 4. pH modulation (Polyelectrolyte chain expansion/contraction):
        ph_dev = abs(ph - 7.4)
        if polymer_name == "Chitosan":
            # Chitosan precipitates at pH > 6.5 due to amine deprotonation
            if ph > 6.5:
                ph_factor = 1.0 + 1.2 * (ph - 6.5)  # aggregate formation
                yield_stress = yield_stress * (1.0 + 1.5 * (ph - 6.5))
            else:
                ph_factor = 1.0
        else:
            ph_factor = max(0.4, 1.0 - 0.15 * (ph_dev ** 1.6))

        K_final = max(0.01, K_temp * ph_factor)
        flow_n_final = max(0.15, min(1.2, flow_n))
        yield_stress_final = max(0.0, yield_stress * ph_factor)

        # Apparent viscosity at resting / zero shear limit (1 s^-1)
        zero_shear_visc = yield_stress_final + K_final * (1.0 ** (flow_n_final - 1.0))

        # Classification
        if flow_n_final < 0.95:
            classification = "Non-Newtonian: Shear-Thinning (Pseudoplastic)"
            rheology_status = "Excellent for bioprinting extrusion; viscosity drops under shear"
        elif flow_n_final <= 1.05:
            classification = "Newtonian Flow"
            rheology_status = "Constant viscosity regardless of shear rate; prone to post-extrusion droplet sagging"
        else:
            classification = "Non-Newtonian: Shear-Thickening (Dilatant)"
            rheology_status = "Viscosity increases under shear; high risk of nozzle jam"

        return {
            "K": round(K_final, 4),
            "n": round(flow_n_final, 4),
            "yield_stress_pa": round(yield_stress_final, 2),
            "zero_shear_viscosity_pa_s": round(zero_shear_visc, 3),
            "classification": classification,
            "status_comment": rheology_status,
        }

    @staticmethod
    def calculate_crosslinking(
        crosslinker_name,
        conc_mm,
        time_min,
        temp_c,
        polymer_name,
        conc_w_v=None,
        ph=7.4,
    ):
        """
        Calculates crosslinking kinetics, network formation degree (0-100%),
        and resulting hydrogel equilibrium elastic modulus (Young's Modulus E and G').

        Incorporates:
        1. Concentration scaling: E ~ (c / c_default)^1.7 in affine hydrogel networks.
        2. Percolation threshold chi_c: Below gel point, G' ~ 0 (liquid sol state).
        3. pH modulation for specific crosslinking mechanisms (e.g. Ca2+/alginate protonation, mTG optimum).
        """
        p_info = POLYMERS_DB.get(polymer_name, POLYMERS_DB["Alginate"])
        c_info = CROSSLINKERS_DB.get(crosslinker_name, CROSSLINKERS_DB["CaCl2"])
        max_mod = c_info["max_modulus_kpa"]
        chi_c = p_info.get("chi_c", 0.12)

        # Use default concentration if not passed
        if conc_w_v is None or conc_w_v <= 0:
            conc_w_v = p_info["conc_default"]

        # pH impact on crosslinking mechanism
        if crosslinker_name == "CaCl2":
            # G-block carboxylates have pKa ~ 3.5; protonated at low pH losing Ca2+ binding capacity
            if ph < 5.0:
                ph_cross_factor = 1.0 / (1.0 + 10.0 ** (3.6 - ph))
            else:
                ph_cross_factor = 1.0
            kd = 18.0  # mM
            conc_term = (conc_mm / (conc_mm + kd)) if conc_mm > 0 else 0.0
            time_term = 1.0 - math.exp(-time_min / 6.0)
            degree = min(1.0, conc_term * time_term * ph_cross_factor)

        elif crosslinker_name == "LAP + 405nm":
            # Radical photo-polymerization: chi = 1 - exp(-k_photo * [LAP] * t)
            t_sec = time_min * 60.0
            ph_cross_factor = max(0.7, 1.0 - 0.05 * abs(ph - 7.4))
            degree = (1.0 - math.exp(-0.025 * conc_mm * t_sec)) * ph_cross_factor

        elif crosslinker_name == "Thermal Gelation":
            # Sigmoidal thermal gel transition
            t_gel = 27.5
            degree = 1.0 / (1.0 + math.exp((temp_c - t_gel) / 1.5))

        elif crosslinker_name == "mTG":
            # Enzymatic isopeptide bond formation; pH optimum around 7.0
            ph_cross_factor = math.exp(-0.5 * ((ph - 7.0) / 0.8) ** 2)
            degree = (
                (conc_mm / (conc_mm + 15.0))
                * (1.0 - math.exp(-time_min / 20.0))
                * ph_cross_factor
            )

        elif crosslinker_name == "Genipin":
            # Nucleophilic amine coupling; optimum around pH 7.5
            ph_cross_factor = math.exp(-0.5 * ((ph - 7.5) / 1.0) ** 2)
            degree = (
                (conc_mm / (conc_mm + 4.0))
                * (1.0 - math.exp(-time_min / 35.0))
                * ph_cross_factor
            )

        else:
            degree = 0.5

        degree = max(0.0, min(0.99, degree))
        degree_pct = round(degree * 100.0, 1)

        # Hydrogel concentration scaling factor: E ~ (c / c_default)^1.7
        c_ratio = max(0.1, conc_w_v / p_info["conc_default"])
        conc_scale = c_ratio ** 1.7

        # Percolation threshold (gel point) model:
        # Below chi_c, isolated clusters do not span space -> G' ~ 0 (sol state).
        # Above chi_c, elasticity emerges via power law (degree - chi_c) / (1 - chi_c)
        if degree <= chi_c:
            g_prime_kpa = 0.05
            youngs_modulus_kpa = 0.1
            gel_state = f"Sol State (Below percolation threshold χ_c={int(chi_c*100)}%; ungelled liquid)"
        else:
            eff_chi = (degree - chi_c) / (1.0 - chi_c)
            g_prime_kpa = round(max_mod * (eff_chi ** 1.8) * conc_scale, 2)
            youngs_modulus_kpa = round(3.0 * g_prime_kpa, 1)

            if degree_pct < 40.0:
                gel_state = "Weakly Crosslinked Hydrogel (Soft matrix; susceptible to mechanical deformation)"
            elif degree_pct <= 85.0:
                gel_state = "Optimal Crosslinked Hydrogel (Balanced elasticity, permeability, and cytocompatibility)"
            else:
                gel_state = "Dense Over-crosslinked Network (Stiff matrix; may impede cell spreading and nutrient diffusion)"

        return {
            "crosslinking_degree_pct": degree_pct,
            "percolation_threshold_pct": round(chi_c * 100.0, 1),
            "storage_modulus_g_prime_kpa": g_prime_kpa,
            "youngs_modulus_kpa": youngs_modulus_kpa,
            "gel_state": gel_state,
            "concentration_scaling_factor": round(conc_scale, 2),
        }

    @staticmethod
    def calculate_extrusion_mechanics(
        nozzle_gauge,
        nozzle_type,
        pressure_kpa,
        syringe_vol_ml,
        K,
        n,
        yield_stress_pa=0.0,
    ):
        """
        Solves the Rabinowitsch-Mooney / Herschel-Bulkley non-Newtonian capillary flow equation
        with physically accurate Bagley entrance loss correction.

        Accounts for:
        1. Bagley entrance loss: Delta_P_total = Delta_P_capillary + Delta_P_entrance.
           delta_P_cap = delta_P_total / (1 + e / (L/D)), where e is Bagley end correction.
        2. Yield stress cutoff: If tau_w <= yield_stress, NO FLOW occurs (Q = 0).
        3. Herschel-Bulkley flow integration when tau_w > yield_stress.
        4. Syringe barrel geometry & dispensing duration from syringe_vol_ml.
        """
        n_info = NOZZLES_DB.get(nozzle_gauge, NOZZLES_DB["21G"])
        D_i_mm = n_info["inner_diam_mm"]
        R_m = (D_i_mm / 2.0) * 1e-3  # radius in meters

        # Capillary length L and Bagley coefficient e
        if nozzle_type == "Conical (Tapered)":
            L_eff_m = 0.007  # 7 mm tapered constricted zone
            bagley_e = 0.5   # Conical tips have minimal contraction vortex losses
        else:
            L_eff_m = 0.0127  # Standard 0.5 inch (12.7 mm) cylindrical blunt needle
            bagley_e = 1.8   # Typical Bagley coefficient (1.5 - 2.5) for blunt needles

        L_over_D = L_eff_m / (2.0 * R_m)

        # Net pressure drop across needle capillary using Bagley correction
        delta_P_total_pa = max(100.0, pressure_kpa * 1e3)
        delta_P_cap_pa = delta_P_total_pa / (1.0 + bagley_e / L_over_D)
        entrance_loss_pa = delta_P_total_pa - delta_P_cap_pa

        # Wall shear stress: tau_w = (Delta_P_cap * R) / (2 * L)
        tau_w = (delta_P_cap_pa * R_m) / (2.0 * L_eff_m)

        # Syringe barrel kinematics (using standard BD Luer-Lok internal diameters)
        barrel_diams = {1.0: 4.78, 3.0: 8.66, 5.0: 12.06, 10.0: 14.50, 20.0: 19.05}
        vol_key = min(barrel_diams.keys(), key=lambda k: abs(k - syringe_vol_ml))
        d_barrel_m = barrel_diams[vol_key] * 1e-3
        a_barrel_m2 = math.pi * ((d_barrel_m / 2.0) ** 2)

        # -------------------------------------------------------------
        # YIELD STRESS CHECK: If tau_w <= yield_stress_pa, NO FLOW
        # -------------------------------------------------------------
        if tau_w <= yield_stress_pa:
            return {
                "inner_diameter_mm": D_i_mm,
                "wall_shear_stress_pa": round(tau_w, 2),
                "wall_shear_stress_kpa": round(tau_w / 1e3, 3),
                "entrance_loss_pa": round(entrance_loss_pa, 1),
                "effective_capillary_pressure_pa": round(delta_P_cap_pa, 1),
                "wall_shear_rate_s_1": 0.0,
                "flow_rate_mm3_s": 0.0,
                "extrusion_velocity_mm_s": 0.0,
                "residence_time_s": 999.0,
                "flow_status": "NO FLOW: Wall shear stress below yield stress (τ_w ≤ τ_y)",
                "flow_active": False,
                "plunger_speed_mm_s": 0.0,
                "dispensing_time_s": 9999.0,
            }

        # -------------------------------------------------------------
        # HERSCHEL-BULKLEY CAPILLARY FLOW INTEGRATION
        # -------------------------------------------------------------
        # Shear rate at wall: gamma_dot_w = ((tau_w - tau_y) / K)^(1/n)
        effective_tau = max(0.01, tau_w - yield_stress_pa)
        gamma_dot_w = (effective_tau / K) ** (1.0 / n)

        # Ratio of yield stress to wall shear stress
        a = yield_stress_pa / tau_w  # in [0, 1)

        # Exact analytical integration for Herschel-Bulkley pipe flow:
        # Q = pi * R^3 * (tau_w/K)^(1/n) * (n/(3n+1)) * (1-a)^((n+1)/n) * [1 + 2n/(2n+1)*a + 2n^2/((n+1)(2n+1))*a^2]
        term1 = math.pi * (R_m ** 3) * ((tau_w / K) ** (1.0 / n)) * (n / (3.0 * n + 1.0))
        term2 = (1.0 - a) ** ((n + 1.0) / n)
        poly = 1.0 + (2.0 * n / (2.0 * n + 1.0)) * a + (2.0 * (n ** 2) / ((n + 1.0) * (2.0 * n + 1.0))) * (a ** 2)
        Q_m3_s = term1 * term2 * poly
        Q_mm3_s = Q_m3_s * 1e9  # mm^3/s (or uL/s)

        # Cross-sectional nozzle area & velocity
        A_m2 = math.pi * (R_m ** 2)
        v_avg_m_s = Q_m3_s / A_m2
        v_avg_mm_s = v_avg_m_s * 1e3

        # Residence time in high-shear zone
        t_res_s = max(0.0005, L_eff_m / max(1e-6, v_avg_m_s))

        # Syringe plunger speed and total dispensing time
        v_plunger_mm_s = (Q_m3_s / a_barrel_m2) * 1e3
        v_dispense_total_mm3 = syringe_vol_ml * 1e3
        t_dispense_s = v_dispense_total_mm3 / max(1e-6, Q_mm3_s)

        # Extrusion rate stability evaluation
        if v_avg_mm_s < 0.5:
            flow_status = "Very Slow / Filament under-extrusion (Risk of tip drying and clogging)"
        elif v_avg_mm_s <= 35.0:
            flow_status = "Ideal Bioprinting Extrusion Window (Stable continuous filament deposition)"
        else:
            flow_status = "High-Speed Jetting / Over-extrusion (Excessive shear stress and filament splashing)"

        return {
            "inner_diameter_mm": D_i_mm,
            "wall_shear_stress_pa": round(tau_w, 2),
            "wall_shear_stress_kpa": round(tau_w / 1e3, 3),
            "entrance_loss_pa": round(entrance_loss_pa, 1),
            "effective_capillary_pressure_pa": round(delta_P_cap_pa, 1),
            "wall_shear_rate_s_1": round(gamma_dot_w, 1),
            "flow_rate_mm3_s": round(Q_mm3_s, 3),
            "extrusion_velocity_mm_s": round(v_avg_mm_s, 2),
            "residence_time_s": round(t_res_s, 4),
            "flow_status": flow_status,
            "flow_active": True,
            "plunger_speed_mm_s": round(v_plunger_mm_s, 3),
            "dispensing_time_s": round(t_dispense_s, 1),
        }

    @staticmethod
    def calculate_cell_safety(
        cell_type,
        cell_density_x10_6,
        wall_shear_stress_pa,
        residence_time_s,
        ph,
        flow_active=True,
    ):
        """
        Calculates cumulative mechanical shear stress dose and models mammalian cell viability.
        Based on Leverett-Nevaril hydrodynamics and Zhang et al. cell membrane rupture models.

        Includes:
        1. Acute immediate mechanical lysis viability (t = 0h).
        2. 72-hour delayed post-extrusion apoptotic trajectory: V(t) = V_0 * exp(-k_apop * t).
        """
        c_info = CELLS_DB.get(cell_type, CELLS_DB["Chondrocytes"])
        tau_crit = c_info["tau_crit_pa"]
        base_k_apop = c_info.get("k_apop", 0.06)

        if not flow_active:
            # Material did not flow; no extrusion shear damage
            return {
                "viability_pct": 97.5,
                "safety_rating": "NO EXTRUSION FLOW",
                "stress_dose_index": 0.0,
                "recommendation": "Bioink did not yield or flow through the nozzle. Increase extrusion pressure or choose a wider nozzle.",
                "disclaimer": "Predictions represent immediate mechanical lysis derived from Leverett-Nevaril hydrodynamic models.",
                "viability_trajectory_72h": [
                    {"hour": h, "viability_pct": 97.5} for h in [0, 6, 12, 24, 48, 72]
                ],
            }

        # Cumulative mechanical stress dose S = (tau_w / tau_crit)^1.15 * (t_res / 1.0)^0.38
        stress_ratio = wall_shear_stress_pa / tau_crit
        dose = (stress_ratio ** 1.15) * (max(1e-4, residence_time_s) ** 0.38)

        # Membrane rupture / damage fraction via Hill function
        damage_fraction = (dose ** 2.2) / (1.0 + (dose ** 2.2))

        # Cell density steric crowding penalty
        if cell_density_x10_6 > 10.0:
            crowd_penalty = min(0.12, 0.02 * (cell_density_x10_6 - 10.0))
        else:
            crowd_penalty = 0.0

        # pH osmotic penalty
        ph_dev = abs(ph - 7.4)
        ph_penalty = min(0.20, 0.05 * (ph_dev ** 2))

        # Baseline unprinted cell viability ~97.5%
        viability = max(
            5.0,
            min(
                98.0,
                97.5 * (1.0 - 0.88 * damage_fraction)
                - (crowd_penalty * 100.0)
                - (ph_penalty * 100.0),
            ),
        )
        viability_pct = round(viability, 1)

        # 72-Hour Post-Print Apoptotic Trajectory
        # Sub-lethally sheared cells exhibit delayed apoptosis over 24-72 hours
        effective_k_apop = base_k_apop * (1.0 + 2.5 * damage_fraction)
        trajectory = []
        for h in [0, 6, 12, 24, 48, 72]:
            t_days = h / 24.0
            v_t = max(2.0, round(viability_pct * math.exp(-effective_k_apop * t_days), 1))
            trajectory.append({"hour": h, "viability_pct": v_t})

        # Qualitative safety rating & clinical recommendations
        if wall_shear_stress_pa < 1500.0 and viability_pct >= 90.0:
            safety_rating = "LOW STRESS / SAFE WINDOW"
            rec = "Current extrusion parameters provide a benign hydrodynamic shear environment with high cellular survival."
        elif wall_shear_stress_pa < 5000.0 and viability_pct >= 80.0:
            safety_rating = "MODERATE STRESS / ACCEPTABLE"
            rec = "Acceptable for robust cell lineages (e.g. fibroblasts, chondrocytes). For sensitive stem cells, consider lowering extrusion pressure by 15-20%."
        elif wall_shear_stress_pa < 9000.0 and viability_pct >= 65.0:
            safety_rating = "HIGH STRESS / CAUTION"
            rec = "Elevated shear causes cell membrane permeabilization and sublethal injury. Switch to a wider nozzle (e.g., 20G/21G) or tapered conical geometry."
        else:
            safety_rating = "SEVERE HYDRODYNAMIC DAMAGE / CRITICAL"
            rec = "Excessive mechanical shear stress exceeds membrane rupture threshold. Immediate nozzle gauge enlargement and pressure reduction mandatory."

        disclaimer = (
            "Viability predictions use a Leverett–Nevaril erythrocyte-derived dose law "
            "extended to nucleated cells with per-lineage τ_crit estimates. Predictions "
            "represent immediate mechanical lysis (0 h) and estimated 72 h apoptotic trajectory."
        )

        return {
            "viability_pct": viability_pct,
            "safety_rating": safety_rating,
            "stress_dose_index": round(dose, 3),
            "recommendation": rec,
            "disclaimer": disclaimer,
            "viability_trajectory_72h": trajectory,
        }

    @staticmethod
    def calculate_scaffold_metrics(
        target_tissue,
        nozzle_gauge,
        layer_height_um,
        strand_spacing_um,
        cell_density_x10_6,
        flow_n,
        yield_stress_pa,
        crosslink_degree_pct,
    ):
        """
        Models 3D scaffold architecture, pore morphology, Ouyang printability index (Pr),
        geometric porosity, and inter-cellular spatial distance.

        Includes both:
        1. Geometric Ouyang circularity index: Pr_ouyang = (pi / 4) * (L_pore^2 / A_pore).
        2. Rheological proxy: Pr_proxy = 0.65 + yield_term + crosslink_term.
        """
        t_info = TISSUES_DB.get(target_tissue, TISSUES_DB["Cartilage"])
        n_info = NOZZLES_DB.get(nozzle_gauge, NOZZLES_DB["21G"])
        d_nozzle_um = n_info["inner_diam_mm"] * 1000.0

        # Die swell / strand spreading factor:
        die_swell = 1.0 + 0.18 * (1.0 - min(0.9, flow_n))
        strand_diam_um = d_nozzle_um * die_swell

        # Layer height and strand spacing in um
        h_layer = max(100.0, layer_height_um)
        s_spacing = max(strand_diam_um + 50.0, strand_spacing_um)

        # Theoretical geometric porosity for orthogonal 0-90° grid scaffold:
        solid_fraction = (math.pi * (strand_diam_um ** 2)) / (4.0 * s_spacing * h_layer)
        porosity_pct = max(15.0, min(92.0, (1.0 - solid_fraction) * 100.0))

        # Pore aperture dimensions
        pore_width_um = max(20.0, s_spacing - strand_diam_um)

        # -------------------------------------------------------------
        # PRINTABILITY INDICES: Geometric Ouyang & Rheological Proxy
        # -------------------------------------------------------------
        # Rheological proxy index:
        yield_term = min(0.20, (yield_stress_pa / 40.0) * 0.20)
        crosslink_term = (crosslink_degree_pct / 100.0) * 0.35
        pr_proxy = round(0.65 + yield_term + crosslink_term, 2)

        # Geometric Ouyang index Pr = (pi/4) * (L_pore^2 / A_pore):
        # Under-gelled bioinks round into corners with radius Rc
        L0 = max(10.0, pore_width_um)
        if pr_proxy < 1.0:
            # Rounding of square corners
            rc = strand_diam_um * 0.5 * (1.0 - pr_proxy)
            rc = min(rc, L0 * 0.45)
            a_pore = (L0 ** 2) - (4.0 - math.pi) * (rc ** 2)
            p_pore = 4.0 * L0 - (8.0 - 2.0 * math.pi) * rc
            pr_geom = round((math.pi / 4.0) * ((p_pore ** 2) / (16.0 * max(1.0, a_pore))), 2)
        elif pr_proxy > 1.05:
            # Over-gelled fracture irregularity
            pr_geom = round(1.0 + (pr_proxy - 1.0) * 0.8, 2)
        else:
            pr_geom = 1.0

        if pr_geom < 0.90:
            pr_assessment = "Under-gelled / Liquid-like (Pore corners round, filament sags, poor shape fidelity)"
        elif pr_geom <= 1.10:
            pr_assessment = "Ideal Printability Window (Sharp orthogonal corners, perfect square pores, high shape fidelity)"
        else:
            pr_assessment = "Over-gelled / Brittle (Irregular fractured strands, lack of inter-layer bonding)"

        # Mean inter-cellular spacing: d_cell = (1 / rho_cells)^(1/3)
        cells_per_um3 = cell_density_x10_6 * 1e-6
        if cells_per_um3 > 0:
            intercell_spacing_um = round((1.0 / cells_per_um3) ** (1.0 / 3.0), 1)
        else:
            intercell_spacing_um = 999.0

        # Tissue target matching diagnostics
        poro_in_range = t_info["target_porosity_min"] <= porosity_pct <= t_info["target_porosity_max"]

        return {
            "strand_diameter_um": round(strand_diam_um, 1),
            "pore_width_um": round(pore_width_um, 1),
            "porosity_pct": round(porosity_pct, 1),
            "printability_pr": pr_geom,
            "printability_proxy": pr_proxy,
            "pr_assessment": pr_assessment,
            "intercell_spacing_um": intercell_spacing_um,
            "target_tissue_porosity_matched": poro_in_range,
        }

    @staticmethod
    def calculate_degradation_and_swelling(polymer_name, crosslink_degree_pct, temp_c, ph):
        """
        Simulates hydrogel mass degradation kinetics and equilibrium swelling across Day 0 to 56.

        Incorporates:
        1. Arrhenius hydrolysis kinetics: k(T) = k_0 * exp[ (Ea/R) * (1/T_ref - 1/T) ] with Ea = 52 kJ/mol.
        2. Polymer-specific pH degradation behavior (e.g. Chitosan accelerated dissolution in acid).
        3. Polymer-specific swelling rate constant k_swell.
        """
        p_info = POLYMERS_DB.get(polymer_name, POLYMERS_DB["Alginate"])
        base_k = p_info["degradation_rate_base"]
        base_swelling = p_info["swelling_equilibrium"]
        k_swell_rate = p_info.get("k_swell", 0.40)

        # Higher crosslinking degree retards degradation and reduces swelling
        cross_factor = max(0.1, 1.0 - 0.75 * (crosslink_degree_pct / 100.0))

        # Arrhenius temperature law: k(T) = k_ref * exp[ (Ea/R) * (1/T_ref - 1/T) ]
        # T_ref = 37 °C = 310.15 K, Ea/R ~ 6250 K (Ea ~ 52 kJ/mol)
        T_k = temp_c + 273.15
        T_ref_k = 37.0 + 273.15
        temp_factor = math.exp(6250.0 * ((1.0 / T_ref_k) - (1.0 / T_k)))

        # Polymer-specific pH degradation factor
        if polymer_name == "Chitosan":
            # Soluble and degrades rapidly at acidic pH; aggregates/precipitates at alkaline pH
            if ph < 6.5:
                ph_factor = 1.0 + 0.4 * (6.5 - ph)
            else:
                ph_factor = max(0.5, 1.0 - 0.2 * (ph - 6.5))
        else:
            # General biopolymers: physiological minimum near 7.4; increases in extreme acid or base
            ph_factor = 1.0 + 0.12 * (abs(ph - 7.4) ** 1.3)

        k_deg = max(1e-4, base_k * cross_factor * temp_factor * ph_factor)
        half_life_days = round(math.log(2.0) / k_deg, 1)

        # Swelling ratio at equilibrium Q_eq (%): Flory-Rehner scaling
        q_equilibrium = round(base_swelling * (1.0 - 0.55 * (crosslink_degree_pct / 100.0)), 1)

        time_points = [0, 7, 14, 21, 28, 56]
        mass_retention = []
        swelling_progression = []

        for t in time_points:
            m_pct = round(100.0 * math.exp(-k_deg * t), 1)
            mass_retention.append({"day": t, "remaining_mass_pct": m_pct})

            s_t = round(q_equilibrium * (1.0 - math.exp(-k_swell_rate * t)), 1) if t > 0 else 0.0
            swelling_progression.append({"day": t, "swelling_ratio_pct": s_t})

        return {
            "degradation_constant_k": round(k_deg, 4),
            "half_life_days": half_life_days,
            "equilibrium_swelling_pct": q_equilibrium,
            "swelling_rate_constant": k_swell_rate,
            "mass_retention_curve": mass_retention,
            "swelling_curve": swelling_progression,
        }

    @staticmethod
    def evaluate_tissue_suitability(
        target_tissue,
        polymer_name,
        cell_type,
        growth_factor,
        rheo_data,
        safety_data,
        scaffold_data,
        crosslink_data,
        deg_data,
    ):
        """
        Comprehensive multi-criteria scoring algorithm across 6 key pillars:
        1. Biochemistry & Cytocompatibility
        2. Rheology & Extrusion Behavior
        3. Printability & Shape Fidelity
        4. Cell Survival Safety
        5. Microstructural Porosity
        6. Degradation & Mechanical Stiffness Matching
        """
        t_info = TISSUES_DB.get(target_tissue, TISSUES_DB["Cartilage"])

        # 1. Biochemistry Score (Cell + GF + Polymer match)
        bio_score = 70.0
        if cell_type in t_info["preferred_cells"]:
            bio_score += 15.0
        if growth_factor in t_info["preferred_gf"]:
            bio_score += 15.0
        if polymer_name in t_info["preferred_polymers"]:
            bio_score += 10.0
        bio_score = min(100.0, bio_score)

        # 2. Rheology Score (Shear thinning n < 0.6 is ideal)
        n = rheo_data["n"]
        if n < 0.55:
            rheo_score = 95.0
        elif n <= 0.75:
            rheo_score = 80.0
        elif n <= 0.95:
            rheo_score = 60.0
        else:
            rheo_score = 30.0

        # 3. Printability Score (Pr proximity to 1.0)
        pr = scaffold_data["printability_pr"]
        pr_dev = abs(pr - 1.0)
        print_score = max(20.0, min(100.0, 100.0 - (pr_dev * 160.0)))

        # 4. Cell Viability Score
        viability = safety_data["viability_pct"]
        cell_score = viability

        # 5. Porosity Score (Proximity to target tissue window)
        porosity = scaffold_data["porosity_pct"]
        if t_info["target_porosity_min"] <= porosity <= t_info["target_porosity_max"]:
            poro_score = 98.0
        else:
            dist = min(
                abs(porosity - t_info["target_porosity_min"]),
                abs(porosity - t_info["target_porosity_max"]),
            )
            poro_score = max(25.0, 95.0 - (dist * 3.5))

        # 6. Degradation & Mechanical Modulus Match
        E = max(0.1, crosslink_data["youngs_modulus_kpa"])
        E_ideal = t_info["target_modulus_ideal"]
        mod_ratio = E / E_ideal if E < E_ideal else E_ideal / E
        mech_score = max(20.0, min(100.0, mod_ratio * 100.0))

        # Degradation match
        hl = deg_data["half_life_days"]
        hl_ideal = t_info["degradation_half_life_ideal_days"]
        hl_ratio = hl / hl_ideal if hl < hl_ideal else hl_ideal / hl
        deg_score = max(25.0, min(100.0, hl_ratio * 100.0))

        mech_deg_combined = (mech_score + deg_score) / 2.0

        # Weighted composite score out of 100
        weights = {
            "biochemistry": 0.20,
            "rheology": 0.15,
            "printability": 0.20,
            "cell_survival": 0.20,
            "porosity": 0.10,
            "mechanics_degradation": 0.15,
        }

        total_score = round(
            bio_score * weights["biochemistry"]
            + rheo_score * weights["rheology"]
            + print_score * weights["printability"]
            + cell_score * weights["cell_survival"]
            + poro_score * weights["porosity"]
            + mech_deg_combined * weights["mechanics_degradation"],
            1,
        )

        if total_score >= 85.0:
            verdict = "EXCELLENT CANDIDATE: Highly recommended for biofabrication trial."
        elif total_score >= 70.0:
            verdict = "SUITABLE FORMULATION: Viable for scaffold printing with minor parameter optimization."
        elif total_score >= 50.0:
            verdict = "BORDERLINE FEASIBILITY: Adjust polymer concentration or nozzle gauge to improve outcome."
        else:
            verdict = "UNFAVORABLE MATCH: Significant incompatibility in rheology, stiffness, or cell viability."

        return {
            "total_score": total_score,
            "verdict": verdict,
            "breakdown": {
                "biochemistry": round(bio_score, 1),
                "rheology": round(rheo_score, 1),
                "printability": round(print_score, 1),
                "cell_survival": round(cell_score, 1),
                "porosity": round(poro_score, 1),
                "mechanics_degradation": round(mech_deg_combined, 1),
            },
        }

    @staticmethod
    def run_inverse_design(
        target_tissue,
        desired_porosity_pct,
        desired_modulus_kpa,
        target_layer_height_um,
    ):
        """
        The Inverse Design AI Engine:
        Searches the multi-dimensional parameter space (polymer, concentration,
        crosslinker concentration, crosslinker time, nozzle gauge, and extrusion pressure)
        to identify the formulation that best meets the desired modulus and porosity.
        """
        t_info = TISSUES_DB.get(target_tissue, TISSUES_DB["Cartilage"])
        candidates = t_info["preferred_polymers"]
        preferred_cell = t_info["preferred_cells"][0]
        preferred_gf = t_info["preferred_gf"][0]

        best_candidate = None
        best_score = -1e9

        test_polymers = candidates if candidates else ["Alginate", "GelMA", "Alginate-Gelatin Blend"]
        test_gauges = ["20G", "21G", "22G", "25G"]
        test_pressures = [60.0, 100.0, 140.0, 180.0, 220.0]

        for poly in test_polymers:
            p_info = POLYMERS_DB[poly]
            c_default = p_info["conc_default"]
            cross_name = p_info["default_crosslinker"]
            c_info = CROSSLINKERS_DB.get(cross_name, CROSSLINKERS_DB["CaCl2"])

            # Concentration search values (0.6x, 0.8x, 1.0x, 1.3x, 1.6x)
            c_vals = [
                round(c_default * f, 1)
                for f in [0.6, 0.8, 1.0, 1.3, 1.6]
                if p_info["conc_min"] <= c_default * f <= p_info["conc_max"]
            ]
            if not c_vals:
                c_vals = [c_default]

            # Crosslinker concentration search values
            c_cross_default = c_info["default_conc_mm"]
            if c_cross_default > 0:
                cross_concs = [
                    round(c_cross_default * 0.5, 1),
                    c_cross_default,
                    round(c_cross_default * 1.5, 1),
                ]
            else:
                cross_concs = [0.0]

            # Crosslinking times
            cross_times = [
                max(0.5, round(c_info["default_time_min"] * 0.5, 1)),
                c_info["default_time_min"],
                round(c_info["default_time_min"] * 1.5, 1),
            ]

            for c_val in c_vals:
                rheo = BioInkEngine.calculate_rheology(
                    poly, c_val, p_info["mw_default_kda"], 25.0, 7.4
                )

                for x_conc in cross_concs:
                    for x_time in cross_times:
                        cross = BioInkEngine.calculate_crosslinking(
                            cross_name, x_conc, x_time, 25.0, poly, conc_w_v=c_val, ph=7.4
                        )

                        for gauge in test_gauges:
                            for p_kpa in test_pressures:
                                ext = BioInkEngine.calculate_extrusion_mechanics(
                                    gauge,
                                    "Conical (Tapered)",
                                    p_kpa,
                                    5.0,
                                    rheo["K"],
                                    rheo["n"],
                                    rheo["yield_stress_pa"],
                                )

                                # Skip configurations with zero flow
                                if not ext["flow_active"] or ext["extrusion_velocity_mm_s"] < 0.2:
                                    continue

                                safety = BioInkEngine.calculate_cell_safety(
                                    preferred_cell,
                                    5.0,
                                    ext["wall_shear_stress_pa"],
                                    ext["residence_time_s"],
                                    7.4,
                                    ext["flow_active"],
                                )

                                scaff = BioInkEngine.calculate_scaffold_metrics(
                                    target_tissue,
                                    gauge,
                                    target_layer_height_um,
                                    1000.0,
                                    5.0,
                                    rheo["n"],
                                    rheo["yield_stress_pa"],
                                    cross["crosslinking_degree_pct"],
                                )

                                deg = BioInkEngine.calculate_degradation_and_swelling(
                                    poly, cross["crosslinking_degree_pct"], 37.0, 7.4
                                )

                                suit = BioInkEngine.evaluate_tissue_suitability(
                                    target_tissue,
                                    poly,
                                    preferred_cell,
                                    preferred_gf,
                                    rheo,
                                    safety,
                                    scaff,
                                    cross,
                                    deg,
                                )

                                # Objective function: maximize suitability while penalizing error to target
                                poro_err = abs(scaff["porosity_pct"] - desired_porosity_pct) / 100.0
                                mod_err = abs(
                                    cross["youngs_modulus_kpa"] - desired_modulus_kpa
                                ) / max(20.0, desired_modulus_kpa)

                                composite_obj = (
                                    suit["total_score"]
                                    - (poro_err * 30.0)
                                    - (mod_err * 30.0)
                                )

                                if composite_obj > best_score:
                                    best_score = composite_obj
                                    best_candidate = {
                                        "polymer": poly,
                                        "concentration_w_v": c_val,
                                        "cell_type": preferred_cell,
                                        "growth_factor": preferred_gf,
                                        "crosslinker": cross_name,
                                        "crosslinker_conc_mm": x_conc,
                                        "crosslinker_time_min": x_time,
                                        "nozzle_gauge": gauge,
                                        "nozzle_type": "Conical (Tapered)",
                                        "extrusion_pressure_kpa": p_kpa,
                                        "layer_height_um": target_layer_height_um,
                                        "strand_spacing_um": round(
                                            scaff["strand_diameter_um"] + (desired_porosity_pct * 10.0), 0
                                        ),
                                        "predicted_score": suit["total_score"],
                                        "predicted_viability_pct": safety["viability_pct"],
                                        "predicted_modulus_kpa": cross["youngs_modulus_kpa"],
                                        "predicted_porosity_pct": scaff["porosity_pct"],
                                        "print_speed_mm_s": round(
                                            max(3.0, min(25.0, ext["extrusion_velocity_mm_s"] * 0.9)), 1
                                        ),
                                        "verdict": suit["verdict"],
                                    }

        return best_candidate

    @staticmethod
    def generate_diagnostics(
        polymer,
        crosslinker,
        rheo,
        ext,
        safety,
        scaff,
        cross,
        target_tissue,
        ph,
    ):
        """
        Generates proactive scientific diagnostics and warnings for anomalous configurations.
        """
        flags = []
        t_info = TISSUES_DB.get(target_tissue, TISSUES_DB["Cartilage"])
        c_info = CROSSLINKERS_DB.get(crosslinker, {})

        # Crosslinker applicability check
        applicable = c_info.get("applicable_polymers", [])
        if applicable and polymer not in applicable:
            flags.append({
                "level": "error",
                "code": "INCOMPATIBLE_CROSSLINKER",
                "title": "Crosslinker Chemistry Incompatibility",
                "message": f"{crosslinker} cannot crosslink {polymer}. Applicable polymers: {', '.join(applicable)}."
            })

        # No flow check
        if not ext.get("flow_active", True) or ext.get("flow_rate_mm3_s", 0) <= 0:
            flags.append({
                "level": "error",
                "code": "NO_FLOW",
                "title": "No Flow Condition (τ_w ≤ τ_y)",
                "message": f"Wall shear stress ({ext['wall_shear_stress_pa']} Pa) is below yield stress ({rheo['yield_stress_pa']} Pa). Increase pressure or nozzle diameter."
            })

        # Shear thickening check
        if rheo["n"] > 1.0:
            flags.append({
                "level": "warning",
                "code": "SHEAR_THICKENING",
                "title": "Dilatant / Shear-Thickening Fluid (n > 1.0)",
                "message": f"Flow behavior index n = {rheo['n']} > 1.0. Viscosity increases with shear rate; high probability of needle clogging."
            })

        # Cell viability warning
        if safety["viability_pct"] < 65.0:
            flags.append({
                "level": "error",
                "code": "LOW_VIABILITY",
                "title": "Severe Cell Lysis (Viability < 65%)",
                "message": f"Predicted viability is {safety['viability_pct']}%. Excessive shear stress ({ext['wall_shear_stress_pa']} Pa) exceeds cell membrane rupture threshold."
            })
        elif safety["viability_pct"] < 80.0:
            flags.append({
                "level": "warning",
                "code": "MODERATE_VIABILITY",
                "title": "Sub-lethal Shear Stress (Viability 65-80%)",
                "message": f"Predicted viability is {safety['viability_pct']}%. Consider conical nozzle or lower pressure."
            })

        # Printability check
        pr = scaff["printability_pr"]
        if pr < 0.85:
            flags.append({
                "level": "warning",
                "code": "UNDER_GELLED",
                "title": "Under-Gelled / Liquid Meniscus (Pr < 0.85)",
                "message": f"Printability index Pr = {pr:.2f}. Filament sags and pores round; poor scaffold shape fidelity."
            })
        elif pr > 1.15:
            flags.append({
                "level": "warning",
                "code": "OVER_GELLED",
                "title": "Over-Gelled / Brittle Filament (Pr > 1.15)",
                "message": f"Printability index Pr = {pr:.2f}. Bioink is over-crosslinked; risk of strand fracture and poor inter-layer bonding."
            })

        # Modulus mismatch
        E = cross["youngs_modulus_kpa"]
        E_min = t_info["target_modulus_kpa_min"]
        E_max = t_info["target_modulus_kpa_max"]
        if E < E_min * 0.33 or E > E_max * 3.0:
            flags.append({
                "level": "warning",
                "code": "STIFFNESS_MISMATCH",
                "title": "Tissue Stiffness Mismatch",
                "message": f"Construct Young's Modulus E = {E} kPa deviates significantly from target tissue requirement ({E_min}-{E_max} kPa for {target_tissue})."
            })

        # pH solubility check
        if polymer == "Chitosan" and ph > 6.5:
            flags.append({
                "level": "warning",
                "code": "CHITOSAN_PRECIPITATION",
                "title": "Chitosan Precipitation Risk (pH > 6.5)",
                "message": f"At pH {ph}, chitosan amines deprotonate causing chain aggregation and precipitation. Use pH ≤ 6.2 for homogeneous solutions."
            })
        elif polymer == "Alginate" and ph < 4.0:
            flags.append({
                "level": "warning",
                "code": "ALGINATE_ACID_GEL",
                "title": "Alginate Acid Precipitation (pH < 4.0)",
                "message": f"At pH {ph}, alginate forms an uncontrolled acid gel with loss of Ca2+ chelation responsiveness."
            })

        return flags

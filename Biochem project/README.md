# BioInkSim: A Biochemical and Rheological Simulator for Polymer-Based Tissue Engineering and 3D Bioprinting

> **Semester Project:** Biochemistry Lab & Biomedical Engineering  
> **Course Scope:** Biochemistry · Tissue Engineering · Rheology · Fluid Dynamics · Biomaterials · Computational Modeling  
> **Application:** Modular Python Package (`bioinksim`), Interactive Web Dashboard (`http://localhost:8501`), Pure Python Calculation Engine, Standalone 3D Viewport  

---

## 1. Project Concept & Executive Summary

In contemporary tissue engineering, 3D bioprinting is far more than mechanical additive manufacturing: it is a complex intersection of **macromolecular biochemistry, non-Newtonian fluid mechanics, cellular mechanobiology, and scaffold micro-architecture**. 

**BioInkSim** is a biophysically grounded simulator that models the complete tissue-engineering cascade:
$$\text{Polymer Selection} \longrightarrow \text{Bioink Formulation} \longrightarrow \text{Rheological Flow} \longrightarrow \text{Crosslinking Kinetics} \longrightarrow \text{Capillary Extrusion} \longrightarrow \text{Hydrodynamic Cell Stress} \longrightarrow \text{3D Scaffold Architecture} \longrightarrow \text{Degradation \& Swelling} \longrightarrow \text{Tissue Suitability Scoring}$$

Rather than presenting static information, **BioInkSim calculates and simulates physical responses in real time** when chemical composition, molecular weight, temperature, pH, cell density, crosslinking agents, and extrusion mechanics change.

```
+---------------------------------------------------------------------------------------------------+
|                                       BIOINKSIM WORKFLOW                                          |
+---------------------------------------------------------------------------------------------------+
|                                                                                                   |
|   [Target Tissue]  --->  [Polymer Hydrogel & Mw]  --->  [Cells & Growth Factors]                  |
|                                                                    |                              |
|                                                                    v                              |
|   [Extrusion Nozzle & Pressure] <--- [Crosslinking Kinetics] <--- [Bioink Formulation (c, T, pH)] |
|               |                                                                                   |
|               v                                                                                   |
|   [Fluid Dynamics (Bagley, Q, tau_w)] ---> [Cell Viability (Nevaril Model + 72h Apoptosis)]       |
|               |                                                                                   |
|               v                                                                                   |
|   [3D Scaffold Architecture (Porosity, Ouyang Pr, Cell Spacing)]                                  |
|               |                                                                                   |
|               v                                                                                   |
|   [Degradation & Swelling (Arrhenius Kinetics)] ---> [Multi-Criteria Tissue Suitability / 100]    |
|               |                                                                                   |
|               v                                                                                   |
|   [RHEOLOGY CALIBRATOR] <---> [UNCERTAINTY QUANTIFICATION] <---> [INVERSE DESIGN AI OPTIMIZER]    |
+---------------------------------------------------------------------------------------------------+
```

---

## 2. Core Physical Models & Governing Equations

### MODULE 1 — Polymer/Bioink Preparation Studio
Formulates the bioink with 7 biopolymers:
* **Sodium Alginate:** Polysaccharide from brown algae; ionic chelation with $\text{Ca}^{2+}$ via cooperative egg-box coordination.
* **Gelatin (Type A):** Denatured collagen; thermo-reversible triple-helix formation ($T_{\text{gel}} \approx 27.5^\circ\text{C}$); abundant native RGD integrin-binding sequences.
* **GelMA (Gelatin Methacryloyl):** Photo-crosslinkable semi-synthetic hydrogel combining thermo-reversible room-temperature extrusion with covalent radical photopolymerization under 405 nm blue / 365 nm UV light (LAP photoinitiator).
* **Collagen Type I:** Mammalian extracellular matrix gold standard; self-assembles into fibrils at $37^\circ\text{C}$ and $\text{pH } 7.4$.
* **Chitosan:** Poly-cationic polysaccharide; antimicrobial and hemostatic; soluble in mild acid ($\text{pH} < 6.2$); aggregates and precipitates at $\text{pH} > 6.5$.
* **Hyaluronic Acid (HA):** Connective tissue glycosaminoglycan (GAG); provides high viscoelastic lubrication and CD44 receptor signaling.
* **Alginate-Gelatin Interpenetrating Network (IPN):** Synergistic composite delivering both high cell attachment (gelatin) and rapid ionic shape fixation (alginate).

---

### MODULE 2 — Non-Newtonian Rheology Simulator
Extrusion bioprinting fundamentally requires **shear-thinning (pseudoplastic)** flow: the bioink must flow smoothly through a narrow gauge needle under shear, then immediately regain high viscosity upon deposition to preserve filament shape fidelity.

#### Governing Equations:
1. **Ostwald–de Waele Power-Law Model:**
   $$\tau = K \dot{\gamma}^n$$
   $$\eta(\dot{\gamma}) = \frac{\tau}{\dot{\gamma}} = K \dot{\gamma}^{n-1}$$
   where:
   * $\tau$ = Shear stress ($\text{Pa}$)
   * $K$ = Consistency index ($\text{Pa}\cdot\text{s}^n$)
   * $\dot{\gamma}$ = Shear rate ($\text{s}^{-1}$)
   * $n$ = Flow behavior index ($n < 1$ shear-thinning; $n = 1$ Newtonian; $n > 1$ dilatant)

2. **Herschel–Bulkley Model (Yield Stress Fluid):**
   $$\tau = \tau_y + K \dot{\gamma}^n \quad (\text{for } \tau > \tau_y)$$
   where $\tau_y$ is the yield stress ($\text{Pa}$). Filaments do not sag or slump post-printing if self-weight stresses remain below $\tau_y$.

3. **Biochemical Modulation of Consistency ($K$) & Flow Index ($n$):**
   * **Concentration Scaling:**
     $$K(c) = K_0 \left(\frac{c}{c_{\text{ref}}}\right)^\alpha \quad (\alpha \approx 2.2 - 3.0)$$
   * **Molecular Weight Dependence:**
     $$K(M_w) = K_c \left(\frac{M_w}{M_{w,\text{ref}}}\right)^{1.4}$$
   * **Temperature Modulation (Arrhenius vs Thermo-Reversible Coil-Helix):**
     * Polysaccharides:
       $$K(T) = K_{\text{ref}} \exp\left[\frac{E_a}{R}\left(\frac{1}{T_K} - \frac{1}{T_{\text{ref},K}}\right)\right]$$
     * Gelatin / GelMA: Sigmoidal physical gel transition around $T_{\text{gel}} = 27.5^\circ\text{C}$.

---

### MODULE 3 — Syringe & Extrusion Hydrodynamics

#### 1. Bagley End Correction for Capillary Flow:
In capillary extrusion through dispensing needles, total pressure drop $\Delta P_{\text{total}}$ consists of capillary viscous dissipation $\Delta P_{\text{cap}}$ and entrance contraction loss $\Delta P_{\text{entrance}}$:
$$\Delta P_{\text{total}} = \Delta P_{\text{cap}} + \Delta P_{\text{entrance}} = 2 \tau_w \left(\frac{L_{\text{eff}}}{R} + e_B \right)$$
where $e_B$ is the Bagley end correction coefficient (typically $e_B \approx 1.5 - 2.5$ in units of $L/D$ for blunt needles, and $e_B \approx 0.4 - 0.7$ for conical tapered tips).

Therefore, effective capillary pressure drop is:
$$\Delta P_{\text{cap}} = \frac{\Delta P_{\text{total}}}{1 + \frac{e_B}{L/D}}$$
For a standard 21G blunt needle ($L/D \approx 24.7$), $\Delta P_{\text{cap}} \approx 0.93 \cdot \Delta P_{\text{total}}$ (entrance loss is ~7%), preventing the 2.4× underestimation of wall shear stress.

Capillary wall shear stress:
$$\tau_w = \frac{\Delta P_{\text{cap}} \cdot R}{2 L_{\text{eff}}}$$

#### 2. Yield Stress Cutoff & Herschel–Bulkley Pipe Flow:
* **No-Flow Condition:** If $\tau_w \le \tau_y$, the material inside the needle remains an unyielded solid plug.
  $$Q = 0, \quad \bar{v} = 0, \quad \dot{\gamma}_w = 0$$
* **Yielded Flow:** When $\tau_w > \tau_y$, let $a = \frac{\tau_y}{\tau_w} \in [0, 1)$. The exact analytical volumetric flow rate for Herschel-Bulkley flow in a cylindrical capillary is:
  $$Q = \pi R^3 \left(\frac{\tau_w}{K}\right)^{1/n} \left(\frac{n}{3n+1}\right) (1-a)^{(n+1)/n} \left[ 1 + \frac{2n}{2n+1} a + \frac{2n^2}{(n+1)(2n+1)} a^2 \right]$$
  * In the limit $\tau_y \to 0$ ($a \to 0$), this simplifies exactly to the Rabinowitsch–Mooney equation.
  * In the Newtonian limit $n = 1, \tau_y = 0$, it reproduces the Hagen–Poiseuille equation: $Q = \frac{\pi R^4 \Delta P_{\text{cap}}}{8 \mu L}$.

#### 3. Syringe Barrel Dynamics:
Using standard BD Luer-Lok barrel diameters ($D_{\text{barrel}} = 12.06\text{ mm}$ for 5 mL syringe), plunger linear speed and dispense run times are calculated:
$$v_{\text{plunger}} = \frac{Q}{\pi (D_{\text{barrel}} / 2)^2}, \quad t_{\text{dispense}} = \frac{V_{\text{syringe}}}{Q}$$

---

### MODULE 4 — Hydrodynamic Shear Stress & Cell Viability

Extrusion exposes encapsulated cells to severe hydrodynamic shear stress at the capillary wall.

#### Biophysical Viability Model:
1. **Leverett–Nevaril Mechanical Stress Dose ($S_{\text{dose}}$):**
   $$S_{\text{dose}} = \left(\frac{\tau_w}{\tau_{\text{crit}}}\right)^{1.15} \cdot \left(\frac{t_{\text{res}}}{1.0\text{ s}}\right)^{0.38}$$
   where $\tau_{\text{crit}}$ is the lineage-specific membrane rupture threshold:
   * Dermal Fibroblasts: $\tau_{\text{crit}} \approx 8000\text{ Pa}$ *(Nevaril et al. 1968, Leverett et al. 1972)*
   * Primary Osteoblasts: $\tau_{\text{crit}} \approx 7200\text{ Pa}$ *(Zhang et al. 2017)*
   * Articular Chondrocytes: $\tau_{\text{crit}} \approx 6200\text{ Pa}$ *(Blaeser et al. 2016)*
   * Endothelial Cells (HUVECs): $\tau_{\text{crit}} \approx 4400\text{ Pa}$ *(Gao et al. 2018)*
   * Mesenchymal Stem Cells (hMSCs): $\tau_{\text{crit}} \approx 3400\text{ Pa}$ *(Ning et al. 2018)*
2. **Acute Mechanical Membrane Rupture Fraction ($D_{\text{fract}}$):**
   $$D_{\text{fract}} = \frac{S_{\text{dose}}^{2.2}}{1.0 + S_{\text{dose}}^{2.2}}$$
3. **Immediate Post-Print Viability ($t = 0\text{ h}$):**
   $$\text{Viability}_{0\text{h}}(\%) = 97.5\% \cdot (1 - 0.88 \cdot D_{\text{fract}}) - \text{CrowdPenalty} - \text{pHPenalty}$$
4. **72-Hour Post-Print Apoptotic Trajectory:**
   Sub-lethal hydrodynamic shear causes delayed apoptotic loss over 24–72 hours:
   $$V(t) = V_{0\text{h}} \cdot \exp\left[ -k_{\text{apop}} \left(1 + 2.5 D_{\text{fract}}\right) \cdot \frac{t}{72\text{ h}} \right]$$

---

### MODULE 5 — 3D Scaffold Architecture & Microstructural Modeling
1. **Die Swell & Strand Diameter ($d_s$):**
   $$d_s = D_i \cdot \left[1 + 0.18 \cdot (1 - \min(0.9, n))\right]$$
2. **Geometric Ouyang Printability Index ($Pr$):**
   Defined by Ouyang et al. (*Biofabrication*, 2016) based on square pore circularity:
   $$Pr = \frac{\pi}{4} \cdot \frac{L_{\text{pore}}^2}{A_{\text{pore}}}$$
   * $Pr = 1.00$: **Ideal square pore** with sharp perpendicular intersection.
   * $Pr < 0.90$: **Under-gelled / Liquid-like** (filaments merge, pore corners round).
   * $Pr > 1.10$: **Over-gelled / Brittle** (strand fracture, poor inter-layer bonding).
3. **Theoretical Geometric Porosity ($P_{\text{geom}}$):**
   $$P_{\text{geom}} = \left[1 - \frac{\pi \cdot d_s^2}{4 \cdot s_{\text{spacing}} \cdot h_{\text{layer}}}\right] \times 100\%$$

---

### MODULE 6 — Crosslinking, Concentration Scaling & Percolation Gel Point
1. **Concentration Scaling in Rubber Elasticity:**
   In affine network theory, hydrogel elastic modulus scales with polymer volume fraction/concentration:
   $$E \propto \left(\frac{c}{c_0}\right)^{1.7}$$
2. **Percolation Threshold (Gel Point $\chi_c$):**
   Real hydrogel gelation requires connectivity percolation (Flory-Stockmayer theory):
   * If $\chi \le \chi_c$: network has not percolated $\implies G' \approx 0.05\text{ kPa}, E \approx 0.1\text{ kPa}$ (sol state).
   * If $\chi > \chi_c$: elasticity emerges as:
     $$G' = G'_{\text{max}} \cdot \left(\frac{\chi - \chi_c}{1 - \chi_c}\right)^{1.8} \cdot \left(\frac{c}{c_0}\right)^{1.7}$$
     $$E \approx 3 \cdot G'$$

---

### MODULE 7 — Polymer Degradation & Swelling Kinetics
1. **Arrhenius Hydrolysis Law:**
   $$k_{\text{deg}}(T) = k_{\text{ref}} \cdot \exp\left[ \frac{E_{a,\text{deg}}}{R} \left(\frac{1}{T_{\text{ref}}} - \frac{1}{T_K}\right) \right]$$
   where $T_{\text{ref}} = 37^\circ\text{C} = 310.15\text{ K}$ and $E_{a,\text{deg}} \approx 52\text{ kJ/mol}$.
2. **Swelling Kinetics with Polymer-Specific Rate Constants:**
   $$Q_s(t) = Q_{\text{eq}} \cdot \left[ 1 - \exp(-k_{\text{swell}} \cdot t) \right]$$
   where $k_{\text{swell}}$ varies by polymer chemistry ($0.85\text{ day}^{-1}$ for Gelatin down to $0.28\text{ day}^{-1}$ for Chitosan).

---

### MODULE 8 — Multi-Criteria Tissue Suitability Matrix

$$\text{Composite Score} = 0.20 S_{\text{bio}} + 0.15 S_{\text{rheo}} + 0.20 S_{\text{print}} + 0.20 S_{\text{cell}} + 0.10 S_{\text{poro}} + 0.15 S_{\text{mech+deg}}$$

---

### MODULE 9 — Rheology Calibration Studio ("Fit to My Data")
Accepts experimental $(\dot{\gamma}, \tau)$ rheometer data and performs pure-Python non-linear regression:
* Fits Ostwald–de Waele: $\tau = K \dot{\gamma}^n$
* Fits Herschel–Bulkley: $\tau = \tau_y + K \dot{\gamma}^n$
* Calculates $R^2$, RMSE, and standard errors of parameters.
* Provides one-click "Apply Calibrated Parameters to Simulator".

---

### MODULE 10 — Uncertainty Quantification (UQ) & Sensitivity Tornado
* **Monte Carlo Propagation:** 500 stochastic draws sampling $K \pm 15\%$, $n \pm 5\%$, $\tau_{\text{crit}} \pm 15\%$, and nozzle diameter $\pm 3\%$. Produces $95\%$ confidence intervals.
* **One-at-a-Time (OAT) Sensitivity Analysis:** Scans parameters $\pm 20\%$ to generate interactive Tornado plots ranking what physical levers move the composite score most.

---

## 3. Project Directory Structure

```
Biochem project/
│
├── bioinksim/                      # Core modular Python package
│   ├── __init__.py                 # Package exports and version metadata
│   ├── data.py                     # POLYMERS_DB, CELLS_DB, BENCHMARKS_DB, etc.
│   ├── engine.py                   # BioInkEngine (physics, rheology, extrusion, UQ)
│   ├── calibration.py              # Pure-Python rheology curve fitting (K, n, tau_y)
│   ├── uq.py                       # Monte Carlo uncertainty quantification engine
│   ├── sensitivity.py              # OAT sensitivity analyzer (Tornado data)
│   ├── server.py                   # Multi-threaded HTTP server with security headers
│   └── frontend/                   # Standalone single-page web dashboard
│       ├── index.html              # HTML5 responsive semantic dashboard
│       ├── style.css               # Modern dark-mode bioprinting stylesheet
│       └── app.js                  # Debounced controller, 3D viewport, charts
│
├── tests/                          # Automated scientific verification test suite
│   ├── __init__.py
│   ├── test_rheology.py            # Shear-thinning, Arrhenius, chitosan pH
│   ├── test_crosslinking.py        # Concentration scaling E~c^1.7, gel point chi_c
│   ├── test_extrusion.py           # Bagley correction, yield stress cutoff
│   ├── test_sanity.py              # Hagen-Poiseuille limit, sol state limits
│   └── test_calibration.py         # Curve fitting accuracy on known data
│
├── cli.py                          # CLI runner, verification runner, benchmarks
├── Bioinksim.py                    # Root backward-compatible entrypoint
└── README.md                       # Comprehensive documentation & equations
```

---

## 4. Software Usage & CLI Commands

### 1. Launch Interactive Web Dashboard
```powershell
python Bioinksim.py
```
* Starts localhost HTTP server at `http://localhost:8501` and opens your default browser.

### 2. Run Automated Verification Test Suite
```powershell
python Bioinksim.py --test
```
* Runs all 10 physical verification modules and closed-form sanity checks.

### 3. Run Standalone CLI Evaluation
```powershell
python Bioinksim.py --cli
```
* Generates an instant command-line report for the selected target tissue and formulation.

### 4. Run Literature Benchmark Parity Report
```powershell
python Bioinksim.py --benchmark
```
* Compares simulator predictions against published literature data (Ouyang 2016, Gao 2018, Müller 2015).

### 5. Download / Export Official PDF Laboratory Report
* **Via Web Dashboard:** Click the **📥 Download PDF Report** button in the header of `http://localhost:8501` to instantly download a publication-grade, 2-page PDF document summarizing the formulation, Bagley hydrodynamics, cellular mechanobiology, 72-hour apoptosis trajectory, and validation audits.
* **Via CLI Command:**
```powershell
python Bioinksim.py --pdf "BioInkSim_Laboratory_Report.pdf"
```
* Generates a standard PDF-1.4 report directly in pure Python with zero external pip dependencies.

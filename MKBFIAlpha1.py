'''
Welcome to MKBFIAlpha1, short for Merged Kostas Beer Fermentation Isothermal: Alpha Edition, Version 1
Compared to MKBFIAlpha1_lite, this Python file:
- ±2 K initial temperature is still not considered, because this is an isothermal BF model
- Adds the shaded region of the EDS in both the graphs (ethanol yield and sugar consumption rates)

In terms of faults:
- The One-At-a-Time Response Curves do not have units, because the code bundles the graphs together.
So, individual y-axes cannot be made in this format. (If fixed, they will be in MKBFIAlpha2)
'''

import math
import numpy as np
import matplotlib.pyplot as plt

# G = glucose, M = maltose, N = maltotriose, X = biomass, E = ethanol, VDK = vicinal diketones, EA = ethyl acetate

# =============
# Parameters
# =============

# MWs [kg/mol]
MW_G = 180.16 / 1000
MW_M = 342.30 / 1000
MW_N = 504.44 / 1000
MW_E = 46.07 / 1000
MW_CO2 = 44.01 / 1000
MW_H2O = 18.015 / 1000

MW_C = 12.01 / 1000
MW_H = 1.008 / 1000
MW_O = 16.00 / 1000
MW_Ni = 14.01 / 1000

# Sugars
G0 = 153.333    # Initial Glucose, mol/m3
M0 = 242.105    # Initial Maltose, mol/m3
N0 = 54.762     # Initial Maltotriose, mol/m3

# convert to kg/m^3
G0_kg_m3 = G0 * MW_G
M0_kg_m3 = M0 * MW_M
N0_kg_m3 = N0 * MW_N

'''
13.1 Plato corresponds to specific gravity of 1.053
~= 138 kg/m3 of total dissolved sugars 
pg. 112 of Stewart, G. G. (2016). Brewing and Distilling Yeasts states:
Sugar       | Percentage Composition
Glucose     | 10-15%
Maltose     | 50-60%
Maltotriose | 15-20%
We're using maple syrup, mainly comprised of sucrose, which breaks down into glucose and fructose
Thus, we'll assume a higher glucose content: G=20%, M=60%, N=20%
G=0.18 kg/mol | M=0.342 kg/mol | N=0.504 kg/mol
i.e. G0 = 0.2 x 138kg / 0.18kg/mol = 153.2 mol/m3
'''

# Initial Values
X0_old_mL = 4.586 # Initial Biomass, million cells/mL
X0_old_m3 = X0_old_mL * 1000000 # Initial Biomass, million cells/m3
# X0 = 0.35 million cells/mL per degree Plato (13.1) = 0.35 * 13.1

# "From this the dry mass of a single yeast cell is calculated as 4.6 x 10^-8 mg." - https://www.tipbiosystems.com/wp-content/uploads/2023/12/AN102-Yeast-Cell-Count_2019_03_17-1.pdf
dry_mass_per_cell_kg = 4.6e-14 # kg per cell
dry_mass_per_million_cells_kg = dry_mass_per_cell_kg * 1000000 # kg per million cells

# 32 percent dry mass % - 68 water mass % - https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/dry-mass
fraction_dry = 0.32

X0_kg_per_m3 = X0_old_m3 * dry_mass_per_million_cells_kg # kg of biomass/m3 of reactor

print(f"Biomass Initial Concentration: {X0_kg_per_m3:.3f} kg/m³")

# "C:H(1.613):O(0.557):N(0.158)" - https://bionumbers.hms.harvard.edu/bionumber.aspx?id=101801&utm

MW_biomass_kg = (5*MW_C + 8*MW_H + 3*MW_O + 1*MW_Ni) # kg/mol

print(f"Biomass Molecular Weight: {MW_biomass_kg:.3f} kg/mol")

X0 = X0_kg_per_m3 / MW_biomass_kg # mol/m3

print(f"Biomass Initial Concentration: {X0:.3f} mol/m³")

E0 = 0.0        # Initial Ethanol, mol/m3
VDK0 = 0.0      # Initial VDKS, ppm
EA0 = 0.0       # Initial Ethyl Acetate, ppm

# Yields
# biomass yields from sugars - in molar ratio - https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
Y_XG = 0.134 # Biomass yield from Glucose
Y_XM = 0.268 # Biomass yield from Maltose
Y_XN = 0.402 # Biomass yield from Maltotriose

# ethanol yields from sugars - in molar ratio - https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
Y_EG = 1.92    # Ethanol yield from Glucose # should be for perfect mass balance 1.777
Y_EM = 3.84     # Ethanol yield from Maltose # should be for perfect mass balance 3.553
Y_EN = 5.76    # Ethanol yield from Maltotriose # should be for perfect mass balance 5.330

CO2_eff = 1 # 0.8675 # efficiency for CO2/ethanol relative to theoretical
CO20 = 0.0       # Initial Carbon Dioxide, mol/m3

'''
Ymolar = Ymass x (MW Sugar / MW Ethanol)

Glucose: C6H1206 -> 2C2H5OH + 2CO2
Ymass = 0.512 g ethanol / g glucose
Ymolar = 2 mol ethanol / mol glucose

Maltose: C12H22O11 + H2O -> 4C2H5OH + 4CO2
Ymass = 0.538 g ethanol / g maltose
Ymolar = 4 mol ethanol / mol maltose

Maltotriose: C18H32016 + 3H2O -> 6C2H5OH + 6CO2
Ymass =  0.548 g ethanol / g maltotriose
Ymolar = 6 mol ethanol / mol maltotriose

This would assume 100% conversion efficiency, generally:
"Among the yeast strains tested, S. cerevisiae UFLA CA 155 had the lowest Ef (77.5%)
S. cerevisiae VR-1 achieved the highest value of Ef (96%) at the end of fermentation"
Citation:
Duarte, W. F.; Dragone, Giuliano; Dias, D. R.; Oliveira, J. M.; Teixeira, J. A.; Schawn, R. F., Efficiency of sugar-to-ethanol conversion by different Saccharomyces cerevisiae strains during raspberry must fermentation. SBFC 2010 - 32nd Symposium on Biotechnology for Fuels and Chemicals. No. 11-23, Clearwater Beach, Florida, USA, April 19-22, 2010
Thus, we'll assume the average
(0.775 + 0.96) / 2 = 0.8675 (87%) * theoretical maximum yields

OR from Gee and Ramiex they are:
YG = 1.92; YM = 3.84; YN = 5.76 [mol/mol]
Which is 96% of theoretical maximum yields
NEEDS DISCUSSED
'''

"""
Vi = Vi0*exp[-Evi/R(T + 273.15)]
Ki = ki0*exp[-Eki/R(T + 273.15)]
Ki' = ki'0*exp[-Eki/R(T + 273.15)]
Reference: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308 

P(T)= P0e[-E/R(T + 273.15)]]
"""

T_ref = 293.15 # K

# simulation conditions
T = 293.15      # Simulation Temperature, K
V = 54.35      # Volume of Liquid, m^3

# Growth Rates
# Arrhenius frequency factor for maximum velocity in h^-1 - Reference: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
V_G0 = np.exp(35.77)
V_M0 = np.exp(16.40)
V_N0 = np.exp(10.59)

# Arrhenius activation energy for maximum velocity in cal/mol - Reference: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
E_VG = 22630
E_VM = 11310
E_VN = 7160

# maximum specific growth rates in h^-1
V_G = V_G0 * np.exp(-E_VG/(1.987*T))  # R = 1.987 cal/mol*K
V_M = V_M0 * np.exp(-E_VM/(1.987*T))
V_N = V_N0 * np.exp(-E_VN/(1.987*T))

#V_G = 0.44     #https://www.sciencedirect.com/science/article/abs/pii/S0141022996001147
#V_M = 0.34     #https://pmc.ncbi.nlm.nih.gov/articles/PMC1214619/
#V_N =

# Michaelis constants in mol/m^3 - Reference: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
K_G = np.exp(-121.3) * np.exp(-(-68600)/(1.987*T))
K_M = np.exp(-19.15) * np.exp(-(-14400)/(1.987*T))
K_N = np.exp(-26.78) * np.exp(-(-19900)/(1.987*T))

# inhibition constants in mol/m^3 - Reference: https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
K_G_prime = np.exp(23.33) * np.exp(-10220/(1.987*T))
K_M_prime = np.exp(55.61) * np.exp(-26350/(1.987*T))

# yeast constants - ignored - model simplified
# X_max = 0     # Max Biomass Density, million cells/mL
# k_d = 0      # Death Rate, h-1

# Vicinal Diketone - https://www.mdpi.com/2304-8158/11/22/3602
Y_VDKX = 0.651 # mass ratio [g_VDK/g_X]
c_VDK = 2.50e-2 # h-1
d_VDK = -5.60e-2 # h-1
k_VDK = c_VDK * math.log(T) + d_VDK # h-1

# Ethyl Acetate (assuming linear extrapolation)
Y_EAS = 1.237e-3 + (T - 14.5) * ((1.237e-3 - 9.870e-4)/(14.5 - 12)) # molar ratio - https://onlinelibrary.wiley.com/doi/epdf/10.1002/j.2050-0416.1994.tb00830.x
Y_EAX = Y_EAS / ((Y_XG + Y_XM + Y_XN)/3) # mole of ethyl acetate formed per mole of biomass consumed

# simulation settings in h
t0 = 0.0        # Start time
tf = 360.0      # End time (15 days)

# solver settings
n_points = 10000
# rel_tol = 1e-8
abs_tol = 1e-10
# max_step = 0.01 # h

"""
FERMENTATION MODEL EQUATIONS - from: https://www.mdpi.com/2304-8158/11/22/3602

1) specific rates:

Equation 1: GLUCOSE SPECIFIC RATE: μ₁ = (V_G × G) / (K_G + G)

Equation 2: MALTOSE SPECIFIC RATE: μ₂ = (V_M × M) / (K_M + M) × K_G′/(K_G′ + G)

Equation 3: MALTOTRIOSE SPECIFIC RATE: μ₃ = (V_N × N) / (K_N + N) × K_G′/(K_G′ + G) × K_M′/(K_M′ + M)

Equation 11: μₓ = Y_XG × μ₁ + Y_XM × μ₂ + Y_XN × μ₃

2) ODES: 

Equation 4: GLUCOSE BALANCE: dG/dt = -μ₁ × X

Equation 5: MALTOSE BALANCE: dM/dt = -μ₂ × X

Equation 6: MALTOTRIOSE BALANCE: dN/dt = -μ₃ × X

Equation 7: BIOMASS BALANCE: dX/dt = (Y_XG × μ1 + Y_XM × μ2 + Y_XN × μ3) × X

Equation 8: ETHANOL BALANCE: dE/dt = Y_EG × (-dG/dt) + Y_EM × (-dM/dt) + Y_EN × (-dN/dt)

Equation 9: VDKs BALANCE: d[VDK]/dt = Y_VDKX × μₓ × X  - k_VDK × [VDK] × X

Equation 10: ETHYL ACETATE BALANCE: d[EA]/dt = Y_EAX × (μ₁ + μ₂ + μ₃) × X
"""


# =============
# Algebraics
# =============

def specific_rates(G, M, N):
    # equation 1
    if G < abs_tol:
        μ1 = 0
    else:
        μ1 = (V_G * G) / (K_G + G)  # h^-1

    # equation 2
    if M < abs_tol:
        μ2 = 0
    else:
        μ2 = (V_M * M) / (K_M + M) * K_G_prime / (K_G_prime + G)  # h^-1

    # equation 3
    if N < abs_tol:
        μ3 = 0
    else:
        μ3 = (V_N * N) / (K_N + N) * K_G_prime / (K_G_prime + G) * K_M_prime / (K_M_prime + M)  # h^-1

    return μ1, μ2, μ3


# =============
# ODEs
# =============

def fermentation_odes(t, y):

    G, M, N, X, E, CO2, VDK, EA = y

    # prevent non-negative concentrations
    G = max(0, G)
    M = max(0, M)
    N = max(N, 0)
    X = max(0, X)

    μ1, μ2, μ3 = specific_rates(G, M, N)

    # equation 4
    if G > abs_tol:
        dG_dt = -μ1 * X  # mol/(m^3*h)
    else:
        dG_dt = 0

    # equation 5
    if M > abs_tol:
        dM_dt = -μ2 * X  # mol/(m^3*h)
    else:
        dM_dt = 0

    # equation 6
    if N > abs_tol:
        dN_dt = -μ3 * X  # mol/(m^3*h)
    else:
        dN_dt = 0

    # equation 7
    # growth_rate = (μ1 + μ2 + μ3) * X * (1 - X / X_max) # million cells/(mL*h)
    # death_rate = k_d * X
    # dX_dt = growth_rate - death_rate
    dX_dt = (Y_XG * μ1 + Y_XM * μ2 + Y_XN * μ3) * X # mol/(m^3*h)

    # equation 8
    dE_dt = Y_EG * (-dG_dt) + Y_EM * (-dM_dt) + Y_EN * (-dN_dt) # mol/(m^3*h)

    # CO2 produced from sugar consumption (TRACKING)
    dCO2_dt = Y_EG * (-dG_dt) + Y_EM * (-dM_dt) + Y_EN * (-dN_dt) # mol/(m^3*h)

    # equation 11
    μx = Y_XG * μ1 + Y_XM * μ2 + Y_XN * μ3  # h^-1

    # equation 9
    dVDK_dt = Y_VDKX * μx * X - k_VDK * VDK * X  # mol/(m^3*h)

    # equation 10
    dEA_dt = Y_EAX * (μ1 + μ2 + μ3) * X  # mol/(m^3*h)

    return np.array([dG_dt, dM_dt, dN_dt, dX_dt, dE_dt, dCO2_dt, dVDK_dt, dEA_dt])


# =============
# RK4 Solver
# =============

def rk4_step(f, t, y, dt):

    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt * k1 / 2)
    k3 = f(t + dt / 2, y + dt * k2 / 2)
    k4 = f(t + dt, y + dt * k3)
    y_new = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    y_new = np.maximum(y_new, 0)
    return y_new

# ===============
# Solve ODEs Helper
# ===============

# time array
t = np.linspace(t0, tf, n_points)
dt = t[1] - t[0]

def solve_odes(G0, M0, N0, X0, E0, CO20, VDK0, EA0):

    solution = np.zeros((n_points, 8))
    solution[0] = [G0, M0, N0, X0, E0, CO20, VDK0, EA0]

    for i in range(n_points - 1):
        solution[i + 1] = rk4_step(
            fermentation_odes, t[i], solution[i], dt
        )

    return t, solution


# =============
# Solve ODEs
# =============

# solution array
t, solution = solve_odes(G0, M0, N0, X0, E0, CO20, VDK0, EA0)

for i in range(n_points - 1):
    solution[i + 1] = rk4_step(fermentation_odes, t[i], solution[i], dt)

glucose = solution[:, 0]
maltose = solution[:, 1]
maltotriose = solution[:, 2]
biomass = solution[:, 3]
ethanol = solution[:, 4]
co2 = solution[:, 5]
vdks = solution[:, 6]
ethyl_acetate = solution[:, 7]

glucose_rate = np.zeros(n_points)
maltose_rate = np.zeros(n_points)
maltotriose_rate = np.zeros(n_points)

def compute_sugar_rates(solution):

    for i in range(n_points):
        glucose_rate[i], maltose_rate[i], maltotriose_rate[i] = specific_rates(
            solution[i, 0],  # glucose
            solution[i, 1],  # maltose
            solution[i, 2],  # maltotriose
        )

    return glucose_rate, maltose_rate, maltotriose_rate

# Sugar bounds
t, solution_s_low = solve_odes(0.95*G0, 0.95*M0, 0.95*N0, X0, E0, CO20, VDK0, EA0)
t, solution_s_high = solve_odes(1.05*G0, 1.05*M0, 1.05*N0, X0, E0, CO20, VDK0, EA0)

ethanol_s_low = solution_s_low[:, 4]
ethanol_s_high = solution_s_high[:, 4]

mu1_s_low, mu2_s_low, mu3_s_low = compute_sugar_rates(solution_s_low)
mu1_s_high, mu2_s_high, mu3_s_high = compute_sugar_rates(solution_s_high)

# =============
# Plotting
# =============

time = t
G = glucose; M = maltose; N = maltotriose
X = biomass; E = ethanol
VDK = vdks; EA = ethyl_acetate
mu1 = glucose_rate; mu2 = maltose_rate; mu3 = maltotriose_rate

# Sugar variables
total_sugar = G + M + N
total_sugar_init = total_sugar[0]
sugar_consumed = total_sugar_init - total_sugar

# Calculate ethanol yield while avoiding division by zero
ethanol_yield = np.divide(E, sugar_consumed, out=np.zeros_like(E), where=sugar_consumed!=0)
ethanol_yield_s_high = np.divide(ethanol_s_high, sugar_consumed, out=np.zeros_like(ethanol_s_high), where=sugar_consumed!=0)
ethanol_yield_s_low = np.divide(ethanol_s_low, sugar_consumed, out=np.zeros_like(ethanol_s_low), where=sugar_consumed!=0)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,15))
plt.subplots_adjust(hspace=0.4)

# Ethanol Yield vs Theoretical Maximum
ax1.plot(time, ethanol_yield, 'b', linewidth=2, label='Actual Yield')
ax1.fill_between(time, ethanol_yield_s_high, ethanol_yield_s_low, alpha=0.2, label='Sugar ±5%')
ax1.axhline(y=0.511, color='r', linestyle='--', linewidth=2, label='Theoretical Max')
ax1.text(0,0.52, 'Theoretical Max (0.511)', color='r', fontsize=10, fontweight='bold')
ax1.set_title('Ethanol Yield vs Theoretical Maximum')
ax1.set_xlabel('Time (h)')
ax1.set_ylabel('Ethanol Yield (mol[ethanol]/mol[sugar])')
ax1.legend(loc='lower right')
ax1.grid(True, alpha=0.3)

# Sugar Consumption Rates
ax2.plot(time, mu1, 'b-', label='Glucose Rate (μ₁)')
ax2.fill_between(time, mu1_s_high, mu1_s_low, alpha=0.2, label='Sugar ±5%')
ax2.plot(time, mu2, 'g-', label='Maltose Rate (μ₂)')
ax2.fill_between(time, mu2_s_high, mu2_s_low, alpha=0.2, label='Sugar ±5%')
ax2.plot(time, mu3, 'r-', label='Maltotriose Rate (μ₃)')
ax2.fill_between(time, mu3_s_high, mu3_s_low, alpha=0.2, label='Sugar ±5%')
ax2.set_title('Specific Sugar Consumption Rates')
ax2.set_ylabel('Specific Rate (h⁻¹)')
ax2.set_xlabel('Time (h)')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.show()

def run_fermentation(
    G0_in,
    M0_in,
    N0_in,
    #T0_in,   # placeholder for later non-isothermal work
):
    """
    Runs the isothermal fermentation model with modified initial conditions.
    Returns scalar outputs for sensitivity analysis.
    """

    # --- override initial conditions ---
    y0 = np.array([
        G0_in,
        M0_in,
        N0_in,
        X0,
        E0,
        CO20,
        VDK0,
        EA0,
    ])

    # time array
    t = np.linspace(t0, tf, n_points)
    dt = t[1] - t[0]

    solution = np.zeros((n_points, 8))
    solution[0] = y0

    for i in range(n_points - 1):
        solution[i + 1] = rk4_step(
            fermentation_odes, t[i], solution[i], dt
        )

    G = solution[:, 0]
    M = solution[:, 1]
    N = solution[:, 2]
    X = solution[:, 3]
    E = solution[:, 4]

    # ---- quantities of interest (QoIs) ----
    final_ethanol = E[-1]

    # fermentation time (99% sugar consumed)
    sugar_total_0 = G0_in + M0_in + N0_in
    sugar_total = G + M + N

    idx = np.where(sugar_total < 0.01 * sugar_total_0)[0]
    fermentation_time = t[idx[0]] if len(idx) > 0 else tf

    max_rate = np.max(-np.gradient(G + M + N, t))

    return {
        "final_ethanol": final_ethanol,
        "fermentation_time": fermentation_time,
        "max_rate": max_rate,
    }

"""
CHAPTER 3 OF 3) SENSITIVITY ANALYSIS FOR ISOTHERMAL BEER FERMENTATION
"""

"""
SUBCHAPTER 3.1) One-At-a-Time Response Curves (Validating Linearity)
"""

# ===============
# Baseline Operating Point
# ===============

G0_base = G0
M0_base = M0
N0_base = N0
T0_base = T   # even if unused

baseline = run_fermentation(G0_base, M0_base, N0_base)
S0_total = G0_base + M0_base + N0_base

# ===============
# External Disturbance Scenario
# ===============

# Sugar mass range: ±5%
sugar_factors = np.linspace(0.95, 1.05, 20)

# Temperature range: ±2 °C (isothermal → flat curves expected)
temperature_values = np.linspace(T0_base - 2, T0_base + 2, 20)

sugar_results = {key: [] for key in baseline.keys()}

for f in sugar_factors:
    out = run_fermentation(
        G0_base * f,
        M0_base * f,
        N0_base * f,
    )
    for key in out:
        sugar_results[key].append(out[key])

temp_results = {key: [] for key in baseline.keys()}

for T in temperature_values:
    out = run_fermentation(
        G0_base,
        M0_base,
        N0_base,
    )
    for key in out:
        temp_results[key].append(out[key])

# ===============
# One-At-a-Time Response Curves
# ===============

def plot_oat_curves(x, results, baseline_val, xlabel, title_prefix):
    n = len(results)
    fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(6, 2.8 * n), sharex=True)

    if n == 1:
        axes = [axes]

    for ax, (key, y_vals) in zip(axes, results.items()):
        ax.plot(x, y_vals, marker='o')
        ax.axhline(baseline_val[key], linestyle='--', linewidth=0.8)
        ax.set_ylabel(key.replace("_", " "))
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel(xlabel)
    fig.suptitle(title_prefix, y=0.98)
    plt.tight_layout()
    plt.show()

plot_oat_curves(
    sugar_factors,
    sugar_results,
    baseline,
    xlabel="Relative sugar mass (–)",
    title_prefix="OAT response to initial sugar mass"
)

plot_oat_curves(
    temperature_values,
    temp_results,
    baseline,
    xlabel="Initial temperature (K)",
    title_prefix="OAT response to temperature (isothermal model)"
)

"""
SUBCHAPTER 3.2) Normalized Local Sensitivities
"""

# ===============
# Normalized Sensitivity
# ===============

def normalized_sensitivity(
    perturb_fn,
    base_output,
    base_input,
    delta_frac,
    output_key,
):
    """
    perturb_fn(factor) → model output dict
    """

    plus = perturb_fn(1 + delta_frac)[output_key]
    minus = perturb_fn(1 - delta_frac)[output_key]

    dy_dx = (plus - minus) / (2 * delta_frac * base_input)

    S_norm = dy_dx * (base_input / base_output)

    return S_norm

# ===============
# Initial Mass of Sugar Disturbance
# ===============

delta_sugar = 0.05
S0_total = G0_base + M0_base + N0_base

def sugar_perturb(factor):
    return run_fermentation(
        G0_base * factor,
        M0_base * factor,
        N0_base * factor,
    )

# ===============
# Initial Temperature Disturbance
# ===============

delta_T = 2.0  # K

def temperature_perturb(factor):
    """
    Factor is ignored because the model is isothermal,
    but the function signature is preserved intentionally.
    """
    return run_fermentation(
        G0_base,
        M0_base,
        N0_base,
    )

sensitivities = {}

# ===============
# Print of Normalized Local Sensitivities
# ===============

print("\nNORMALIZED LOCAL SENSITIVITIES\n")

for output_key, base_val in baseline.items():

    S_sugar = normalized_sensitivity(
        sugar_perturb,
        base_val,
        S0_total,
        delta_sugar,
        output_key,
    )

    # temperature sensitivity (expected ~0)
    S_temp = (
        temperature_perturb(1.0)[output_key]
        - temperature_perturb(1.0)[output_key]
    ) / (2 * delta_T)

    # normalized form (safe division)
    S_temp_norm = 0.0 if abs(base_val) < 1e-12 else S_temp * (T0_base / base_val)

    sensitivities[output_key] = {
        "Sugar mass": S_sugar,
        "Temperature": S_temp_norm,
    }

    print(f"{output_key}:")
    print(f"  Sugar mass sensitivity : {S_sugar:+.3f}")
    print(f"  Temperature sensitivity: {S_temp_norm:+.3f}")

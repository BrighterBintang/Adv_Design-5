'''
Welcome to MKBFIFinal, short for Merged Kostas Beer Fermentation Isothermal: Final Edition
Compared to MKBFIRelease1, this Python file:
- Has the completed graphs with both EDS showing up on all graphs.
- Intended for submission of the ACED5 Inception Report

Note:
- Graph 5 (Sugar Consumption Kinetics) is broken, so that graph will have to be generated without sensitivity analysis.
- Graph 1's legend is quite large which blocks the line and shaded areas of maltose. The legend needs to be rescaled.
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

MW_VDK = 4 * MW_C + 6 * MW_H + 2 * MW_O # diacetyl

print(f"Molecular weight of VDK: {MW_VDK:.6f} kg/mol")

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

yield_eff = 0.8675 # 0.8675 # efficiency for CO2/ethanol relative to theoretical

# ethanol yields from sugars - in molar ratio - https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.260310308
# Y_EG = 1.92    # Ethanol yield from Glucose # should be for perfect mass balance 1.777
# Y_EM = 3.84     # Ethanol yield from Maltose # should be for perfect mass balance 3.553
# Y_EN = 5.76    # Ethanol yield from Maltotriose # should be for perfect mass balance 5.330

# theoretical maximum: Y_EG = 2, Y_EM = 4, Y_EN = 6 - current evaluations assume 96 percent efficiency (maximum)

Y_EG = 2.0 * yield_eff
Y_EM = 4.0 * yield_eff
Y_EN = 6.0 * yield_eff

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
Y_VDKX_mass = 0.651 # mass ratio [g_VDK/g_X]
Y_VDKX = Y_VDKX_mass * (MW_biomass_kg / MW_VDK) # molar ratio
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
# Update Kinetics (useful for the Effects of Initial Temperature on Output Variables)
# ===============

def set_temp_kinetics(T_input):
    global V_G, V_M, V_N, K_G, K_M, K_N, K_G_prime, K_M_prime

    V_G = V_G0 * np.exp(-E_VG / (1.987 * T_input))  # R = 1.987 cal/mol*K
    V_M = V_M0 * np.exp(-E_VM / (1.987 * T_input))
    V_N = V_N0 * np.exp(-E_VN / (1.987 * T_input))

    K_G = np.exp(-121.3) * np.exp(-(-68600) / (1.987 * T_input))
    K_M = np.exp(-19.15) * np.exp(-(-14400) / (1.987 * T_input))
    K_N = np.exp(-26.78) * np.exp(-(-19900) / (1.987 * T_input))

    K_G_prime = np.exp(23.33) * np.exp(-10220 / (1.987 * T_input))
    K_M_prime = np.exp(55.61) * np.exp(-26350 / (1.987 * T_input))

# =============
# Solve ODEs
# =============

# time array
t = np.linspace(t0, tf, n_points)
dt = t[1] - t[0]

# solution array
solution = np.zeros((n_points, 8))
solution[0] = [G0, M0, N0, X0, E0, CO20, VDK0, EA0]

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

for i in range(n_points):
    glucose_rate[i], maltose_rate[i], maltotriose_rate[i] = specific_rates(glucose[i], maltose[i], maltotriose[i])

# ===============
# Solve ODEs for External Disturbance Scenario
# ===============

y0_standard = np.array([G0, M0, N0, X0, E0, CO20, VDK0, EA0])
T0_standard = T

def perturb_initial_conditions(
    y0,
    T0,
    sugar_factor=1.0,
    delta_T=0.0
):
    y0_new = y0.copy()

    # sugars: glucose, maltose, maltotriose
    y0_new[0:3] *= sugar_factor

    T0_new = T0 + delta_T

    return y0_new, T0_new

# +5% sugars, 0 K
y0_plus, T0_standard = perturb_initial_conditions(
    y0_standard, T0_standard,
    sugar_factor=1.05,
    delta_T=0.0
)

# −5% sugars, 0 K
y0_minus, T0_standard = perturb_initial_conditions(
    y0_standard, T0_standard,
    sugar_factor=0.95,
    delta_T=0.0
)

# +0% sugars, +2 K
y0_standard, T_plus = perturb_initial_conditions(
    y0_standard, T0_standard,
    sugar_factor=1.0,
    delta_T=2.0
)

# +0% sugars, −2 K
y0_standard, T_minus = perturb_initial_conditions(
    y0_standard, T0_standard,
    sugar_factor=1.0,
    delta_T=-2.0
)

def run_simulation(y0, T0):
    set_temp_kinetics(T0)

    solution = np.zeros((n_points, 8))
    solution[0] = y0

    for i in range(n_points - 1):
        solution[i + 1] = rk4_step(
            fermentation_odes,
            t[i],
            solution[i],
            dt
        )

    return solution

sol_standard = run_simulation(y0_standard, T0_standard)
print(sol_standard)
sol_sugar_plus = run_simulation(y0_plus, T0_standard)
print(sol_sugar_plus)
sol_sugar_minus = run_simulation(y0_minus, T0_standard)
print(sol_sugar_minus)
sol_T_plus = run_simulation(y0_standard, T_plus)
print(sol_T_plus)
sol_T_minus = run_simulation(y0_standard, T_minus)
print(sol_T_minus)

def compute_sugar_rates(solution):

    for i in range(n_points):
        glucose_rate[i], maltose_rate[i], maltotriose_rate[i] = specific_rates(
            solution[i, 0],  # glucose
            solution[i, 1],  # maltose
            solution[i, 2],  # maltotriose
        )

    return glucose_rate, maltose_rate, maltotriose_rate






# =============
# Plotting
# =============

time = t
G, M, N = glucose, maltose, maltotriose
X, E = biomass, ethanol
VDK, EA = vdks, ethyl_acetate
mu1, mu2, mu3 = glucose_rate, maltose_rate, maltotriose_rate

# =============
# Plotting but for External Disturbance Scenario
# =============

# +5% Sugar

G_sp, M_sp, N_sp = sol_sugar_plus[:, 0], sol_sugar_plus[:, 1], sol_sugar_plus[:, 2]
X_sp, E_sp = sol_sugar_plus[:, 3], sol_sugar_plus[:, 4]
VDK_sp, EA_sp = sol_sugar_plus[:, 6], sol_sugar_plus[:, 7]
mu1_sp, mu2_sp, mu3_sp = compute_sugar_rates(sol_sugar_plus)

# -5% Sugar

G_sm, M_sm, N_sm = sol_sugar_minus[:, 0], sol_sugar_minus[:, 1], sol_sugar_minus[:, 2]
X_sm, E_sm = sol_sugar_minus[:, 3], sol_sugar_minus[:, 4]
VDK_sm, EA_sm = sol_sugar_minus[:, 6], sol_sugar_minus[:, 7]
mu1_sm, mu2_sm, mu3_sm = compute_sugar_rates(sol_sugar_minus)

# +2 K Temperature

G_Tp, M_Tp, N_Tp = sol_T_plus[:, 0], sol_T_plus[:, 1], sol_T_plus[:, 2]
X_Tp, E_Tp = sol_T_plus[:, 3], sol_T_plus[:, 4]
VDK_Tp, EA_Tp = sol_T_plus[:, 6], sol_T_plus[:, 7]
mu1_Tp, mu2_Tp, mu3_Tp = compute_sugar_rates(sol_T_plus)

# -2 K Temperature

G_Tm, M_Tm, N_Tm = sol_T_minus[:, 0], sol_T_minus[:, 1], sol_T_minus[:, 2]
X_Tm, E_Tm = sol_T_minus[:, 3], sol_T_minus[:, 4]
VDK_Tm, EA_Tm = sol_T_minus[:, 6], sol_T_minus[:, 7]
mu1_Tm, mu2_Tm, mu3_Tm = compute_sugar_rates(sol_T_minus)

# these are the fonts and styles just leave them
plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['axes.labelsize'] = 24
plt.rcParams['axes.titlesize'] = 28
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['lines.linewidth'] = 4.0
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['legend.framealpha'] = 0.95
plt.rcParams['legend.edgecolor'] = 'black'
plt.rcParams['legend.fancybox'] = True
plt.rcParams['legend.shadow'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

# index for glucose depletion
glucose_threshold = 1
depletion_idx = next((i for i, g in enumerate(G) if g <= glucose_threshold), None)
depletion_time = time[depletion_idx] if depletion_idx is not None else None

# colour palette for each variable - using help from https://www.websitecolorchooser.com/business-color-scheme-explorer/
# leave as "colors" and "color" to not bug the package
colors = {'glucose': '#3B82F6', 'maltose': '#52796F', 'maltotriose': '#E53E3E', 'biomass': '#17202A', 'ethanol': '#FF6F3C', 'co2': '#667EEA', 'vdks': '#FFB3C1', 'ethyl_acetate': '#F4A259', 'depletion_line': '#2D3436'}

# vertical dotted line to indicate glucose depletion (and hide ugly transition)
def plot_depletion_line(ax, y_factor=0.9):
    if depletion_idx is not None:
        ax.axvline(depletion_time, color=colors['depletion_line'], linestyle='--', linewidth=2.5, alpha=0.7)
        ax.text(depletion_time, ax.get_ylim()[1]*y_factor, 'Glucose\nDepleted', ha='center', va='top', fontsize=12, fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))

# 1) Sugar Concentrations vs Time
fig1, ax1 = plt.subplots(figsize=(14, 10))

# plot glucose truncated at depletion
ax1.plot(time[:depletion_idx+1] if depletion_idx else time, G[:depletion_idx+1] if depletion_idx else G, color=colors['glucose'], label='Glucose', linewidth=4)
ax1.fill_between(time, G_sp, G_sm, alpha=0.2, color= (0, 0, 1, 0.2), label='Sugar ±5%')
ax1.fill_between(time, G_Tp, G_Tm, alpha=0.2, color= (0, 0.5, 1, 0.2), label='Temperature ±2 K')

# plot others fully
ax1.plot(time, M, color=colors['maltose'], label='Maltose', linewidth=4)
ax1.fill_between(time, M_sp, M_sm, alpha=0.2, color= (0, 1, 0, 0.2), label='Sugar ±5%')
ax1.fill_between(time, M_Tp, M_Tm, alpha=0.2, color= (0.5, 1, 0.5, 0.2), label='Temperature ±2 K')

ax1.plot(time, N, color=colors['maltotriose'], label='Maltotriose', linewidth=4)
ax1.fill_between(time, N_sp, N_sm, alpha=0.2, color= (1, 0, 0, 0.2), label='Sugar ±5%')
ax1.fill_between(time, N_Tp, N_Tm, alpha=0.2, color= (1, 0.5, 0, 0.2), label='Temperature ±2 K')

# glucose depletion line
plot_depletion_line(ax1)

ax1.set(xlabel='Time (h)', ylabel='Concentration (mol/m³)', title='Sugar Consumption Profile')
ax1.legend(loc='upper right', framealpha=0.95)
ax1.grid(True, which='both', alpha=0.2, linestyle=':')
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)

plt.tight_layout()
plt.show()

# 2) Biomass Concentration vs Time
fig2, ax2 = plt.subplots(figsize=(14, 8))

# convert biomass to million cells/mL (make it more clear)
biomass_million_cells_ml = X * MW_biomass_kg / dry_mass_per_million_cells_kg / 1e6

bmc_sp = X_sp * MW_biomass_kg / dry_mass_per_million_cells_kg / 1e6
bmc_sm = X_sm * MW_biomass_kg / dry_mass_per_million_cells_kg / 1e6
bmc_Tp = X_Tp * MW_biomass_kg / dry_mass_per_million_cells_kg / 1e6
bmc_Tm = X_Tm * MW_biomass_kg / dry_mass_per_million_cells_kg / 1e6

ax2.plot(time, biomass_million_cells_ml, color=colors['biomass'], label='Biomass', linewidth=3.5, alpha=0.9)
ax2.fill_between(time, bmc_sp, bmc_sm, alpha=0.2, color= (0, 0, 0, 0.2), label='Sugar ±5%')
ax2.fill_between(time, bmc_Tp, bmc_Tm, alpha=0.2, color= (0.5, 0.5, 0.5, 0.2), label='Temperature ±2 K')

plot_depletion_line(ax2, y_factor=0.95)

ax2.set(xlabel='Time (h)', ylabel='Biomass (million cells/mL)', title='Biomass Growth Profile')
ax2.legend(loc='upper left', framealpha=0.95)

# create twin axis for raw biomass (one in mol/m3 and one in million cells/mL)
# ax2b = ax2.twinx()
# ax2b.plot(time, X, color=colors['biomass'], linestyle=':', linewidth=2, alpha=0.6)
# ax2b.set_ylabel('Biomass (mol/m³)', color=colors['biomass'])
# ax2b.tick_params(axis='y', labelcolor=colors['biomass'])

ax2.grid(True, which='both', alpha=0.1, linestyle=':')
ax2.minorticks_on()

plt.tight_layout()
plt.show()

# 3) Ethanol and CO2 Concentrations vs Time
fig3, ax3 = plt.subplots(figsize=(14, 8))

ax3.plot(time, E, color=colors['ethanol'], label='Ethanol', linewidth=3.5, alpha=0.9)
ax3.fill_between(time, E_sp, E_sm, alpha=0.2, color= (1, 0.5, 0, 0.2), label='Sugar ±5%')
ax3.fill_between(time, E_Tp, E_Tm, alpha=0.2, color= (1, 0.75, 0.25, 0.2), label='Temperature ±2 K')

ax3.set_xlabel('Time (h)')
ax3.set_ylabel('Ethanol (mol/m³)', color=colors['ethanol'])
ax3.tick_params(axis='y', labelcolor=colors['ethanol'])

# CO2 twin axis
# ax3b = ax3.twinx()
# ax3b.plot(time, co2, color=colors['co2'], linestyle='--', linewidth=2.5, alpha=0.8, label='CO₂')
# ax3b.set_ylabel('CO₂ (mol/m³)', color=colors['co2'])
# ax3b.tick_params(axis='y', labelcolor=colors['co2'])

plot_depletion_line(ax3, y_factor=0.8)

lines, labels = ax3.get_legend_handles_labels()
# lines2, labels2 = ax3b.get_legend_handles_labels()
ax3.legend(lines, labels, loc='upper left', framealpha=0.95)

ax3.set_title('Ethanol Production Profile')
ax3.grid(True, which='both', alpha=0.1, linestyle=':')
ax3.minorticks_on()

plt.tight_layout()
plt.show()

# 4) Flavour Compounds Concentrations vs Time
fig4, ax4 = plt.subplots(figsize=(14, 8))

ax4.plot(time, VDK, color=colors['vdks'], label='VDKs', linewidth=3, alpha=0.9)
ax4.fill_between(time, VDK_sp, VDK_sm, alpha=0.2, color= (1, 0.75, 0.8, 0.2), label='Sugar ±5%')
ax4.fill_between(time, VDK_Tp, VDK_Tm, alpha=0.2, color= (0.8, 0.5, 0.8, 0.2), label='Temperature ±2 K')

ax4.plot(time, EA, color=colors['ethyl_acetate'], label='Ethyl Acetate', linewidth=3, alpha=0.9)
ax4.fill_between(time, EA_sp, EA_sm, alpha=0.2, color= (1, 0.5, 0, 0.2), label='Sugar ±5%')
ax4.fill_between(time, EA_Tp, EA_Tm, alpha=0.2, color= (0.6, 0.3, 0, 0.2), label='Temperature ±2 K')

plot_depletion_line(ax4)

# highlights the values at depletion
ax4.scatter([depletion_time], [VDK[depletion_idx]], color=colors['vdks'], s=100, edgecolor='black', zorder=5)
ax4.scatter([depletion_time], [EA[depletion_idx]], color=colors['ethyl_acetate'], s=100, edgecolor='black', zorder=5)

ax4.set(xlabel='Time (h)', ylabel='Concentration (ppm)', title='Flavour Compounds Production')
ax4.legend(loc='upper left', framealpha=0.95)
ax4.grid(True, which='both', alpha=0.1, linestyle=':')
ax4.minorticks_on()

plt.tight_layout()
plt.show()

# 5) Sugar Consumption Rates vs Time
fig5, ax5 = plt.subplots(figsize=(14, 8))

# glucose rate becomes truncated at depletion
ax5.plot(time[:depletion_idx+1] if depletion_idx else time, mu1[:depletion_idx+1] if depletion_idx else mu1, color=colors['glucose'], label='μ₁ (Glucose)', linewidth=3, alpha=0.9)
ax5.fill_between(time, mu1_sp, mu1_sm, color = (0, 0, 1, 0.2), alpha=0.2, label='Sugar ±5%')
ax5.fill_between(time, mu1_Tp, mu1_Tm, color = (0, 0.5, 1, 0.2), alpha=0.2, label='Temperature ±2 K')

ax5.plot(time, mu2, color=colors['maltose'], label='μ₂ (Maltose)', linewidth=3, alpha=0.9)
ax5.fill_between(time, mu2_sp, mu2_sm, color = (0, 1, 0, 0.2), alpha=0.2, label='Sugar ±5%')
ax5.fill_between(time, mu2_Tp, mu2_Tm, color = (0.5, 1, 0.5, 0.2), alpha=0.2, label='Temperature ±2 K')

ax5.plot(time, mu3, color=colors['maltotriose'], label='μ₃ (Maltotriose)', linewidth=3, alpha=0.9)
ax5.fill_between(time, mu3_sp, mu3_sm, color = (1, 0, 0, 0.2), alpha=0.2, label='Sugar ±5%')
ax5.fill_between(time, mu3_Tp, mu3_Tm, color = (1, 0.5, 0, 0.2), alpha=0.2, label='Temperature ±2 K')

plot_depletion_line(ax5)

# placement for legend is weird so will be manually set not overlap label
ax5.set(xlabel='Time (h)', ylabel='Specific Consumption Rate (h⁻¹)', title='Sugar Consumption Kinetics')
ax5.legend(framealpha=0.95, loc='upper left', bbox_to_anchor=(0.71, 0.8))
ax5.grid(True, which='both', alpha=0.1, linestyle=':')
ax5.minorticks_on()
ax5.set_xlim(left=0)
ax5.set_ylim(bottom=0)

plt.tight_layout()
plt.show()

# 6) Ethanol Yield vs Theoretical Maximum
fig6, ax6 = plt.subplots(figsize=(14, 10))

# calculations for ethanol yield
ethanol_yield = E / (2 * (G0 - G) + 4 * (M0 - M) + 6 * (N0 - N) + 1e-10)

ey_sp = E_sp / (2 * (1.05*G0 - G_sp) + 4 * (1.05*M0 - M_sp) + 6 * (1.05*N0 - N_sp) + 1e-10)
ey_sm = E_sm / (2 * (0.95*G0 - G_sm) + 4 * (0.95*M0 - M_sm) + 6 * (0.95*N0 - N_sm) + 1e-10)
ey_Tp = E_Tp / (2 * (G0 - G_Tp) + 4 * (M0 - M_Tp) + 6 * (N0 - N_Tp) + 1e-10)
ey_Tm = E_Tm / (2 * (G0 - G_Tm) + 4 * (M0 - M_Tm) + 6 * (N0 - N_Tm) + 1e-10)

Y_theoretical_max = 1.0      # 100% theoretical
Y_practical_max = 0.96       # 96% maximum from literature
Y_avg_efficiency = 0.8675    # 86.75% average (what we use)

# plot actual yield
ax6.plot(time, ethanol_yield, 'b-', linewidth=4, label='Actual Yield')
ax6.fill_between(time, ey_sp, ey_sm, alpha=0.2, color= (0, 0, 1, 0.2), label='Sugar ±5%')
ax6.fill_between(time, ey_Tp, ey_Tm, alpha=0.2, color= (0, 1, 1, 0.2), label='Temperature ±2 K')

# plot reference yields
ax6.axhline(Y_theoretical_max, color='r', linestyle='--', linewidth=4, alpha=0.7, label='Theoretical Max (100%)')
ax6.axhline(Y_practical_max, color='g', linestyle='--', linewidth=4, alpha=0.7, label='Practical Max (96%)')
ax6.axhline(Y_avg_efficiency, color='orange', linestyle='--', linewidth=4, alpha=0.7, label='Average Efficiency (86.75%)')

ax6.set(xlabel='Time (h)', ylabel='Ethanol Yield (mol[ethanol]/mol[sugar])', title='Ethanol Yield vs Theoretical and Practical Maximum')
ax6.legend(loc='lower right', framealpha=0.95)
ax6.grid(True, alpha=0.3, linestyle='--')

ax6.set_xlim(left=0, right=time[-1])

ax6.set_ylim(bottom=0.7, top=1.02)

# annotate theoretical max
ax6.text(0.02, 0.95, 'Theoretical Maximum', transform=ax6.transAxes, color='r', fontsize=18, verticalalignment='bottom', fontweight='bold')

# annotate practical max
ax6.text(0.02, 0.83, 'Practical Maximum', transform=ax6.transAxes, color='g', fontsize=18, verticalalignment='bottom', fontweight='bold')

ax6.xaxis.set_major_locator(plt.MaxNLocator(8))
ax6.yaxis.set_major_locator(plt.MaxNLocator(8))

plt.tight_layout()
plt.show()

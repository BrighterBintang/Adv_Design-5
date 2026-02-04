'''
Welcome to MKBFI, short for Merged Kostas Beer Fermentation Isothermal
This Python file includes:
- Kostas's Parameters.py and BFIsothermal.py (formerly Beer Fermentation Isothermal.py [Lines 1-176])
- Bintang's Additions to BFIsothermal.py [Lines 177-233] and STAIsothermal.py
(an adaptation of Sensitivity Analysis Demo.py for Kostas's BFIsothermal code)
'''

import numpy as np
import csv
import matplotlib.pyplot as plt

"""
CHAPTER 1 OF 3) PARAMETERS

PARAMETERS FOR ISOTHERMAL BEER FERMENTATION MODEL

THESE VALUES ARE ALL MADE UP AND NEED TO BE REPLACED WITH REAL VALUES FROM LITERATURE
"""

# G = glucose, M = maltose, N = maltotriose, X = biomass, E = ethanol, VDK = vicinal diketones, EA = ethyl acetate
# initial conditions in mol/m^3
G0 = 100.0
M0 = 100.0
N0 = 100.0
X0 = 0.35       # million cells/mL
E0 = 0.0
VDK0 = 0.0      # ppm
EA0 = 0.0       # ppm

# biomass yields from sugars
Y_XG = 0.1
Y_XM = 0.1
Y_XN = 0.1

# ethanol yields from sugars
Y_EG = 0.5
Y_EM = 0.5
Y_EN = 0.5

# maximum specific growth rates in h^-1
V_G = 0.3
V_M = 0.2
V_N = 0.15

# Michaelis constants in mol/m^3
K_G = 1.0
K_M = 1.0
K_N = 1.0

# inhibition constants in mol/m^3
K_G_prime = 10.0
K_M_prime = 10.0

# yeast constants
X_max = 50.0       # million cells/mL
k_d = 0.01         # h-1

# flavourings constants
k_f = 0.001        # ppm/(h*million cells/mL)
k_r = 0.005        # h-1
k_EA = 0.0001      # ppm/(h*million cells/mL*mol/m3)
k_hyd = 0.0005     # h-1

# simulation conditions
T0 = 293.15      # K
V = 100.0       # m^3

# simulation settings in h
t0 = 0.0
tf = 336.0

# solver settings
n_points = 10000
# rel_tol = 1e-8
# abs_tol = 1e-10
# max_step = 0.01 # h

"""
CHAPTER 2 OF 3) ISOTHERMAL BEER FERMENTATION SOLVER

ISOTHERMAL FERMENTATION SOLVER

Code Logic:
1) Import parameters [already in the same file]
2) Perform algebraic equations 
3) Define ODEs
4) Define RK4
5) Solve ODEs
6) Save to CSV (for plotting)
"""

"""
FERMENTATION MODEL EQUATIONS

1) specific rates:

Equation 1: GLUCOSE SPECIFIC RATE: μ₁ = (V_G × G) / (K_G + G)

Equation 2: MALTOSE SPECIFIC RATE: μ₂ = (V_M × M) / (K_M + M) × K_G′/(K_G′ + G)

Equation 3: MALTOTRIOSE SPECIFIC RATE: μ₃ = (V_N × N) / (K_N + N) × K_G′/(K_G′ + G) × K_M′/(K_M′ + M)

2) ODES: 

Equation 4: GLUCOSE BALANCE: dG/dt = -μ₁ × X

Equation 5: MALTOSE BALANCE: dM/dt = -μ₂ × X

Equation 6: MALTOTRIOSE BALANCE: dN/dt = -μ₃ × X

Equation 7: BIOMASS BALANCE: dX/dt = (μ₁ + μ₂ + μ₃) × X × (1 - X/X_max) - k_d × X

Equation 8: ETHANOL BALANCE: dE/dt = Y_EG × (-dG/dt) + Y_EM × (-dM/dt) + Y_EN × (-dN/dt)

Equation 9: VDKs BALANCE: d[VDK]/dt = k_f × X - k_r × [VDK]

Equation 10: ETHYL ACETATE BALANCE: d[EA]/dt = k_EA × X × E - k_hyd × [EA]
"""


# =============
# Algebraics
# =============

def specific_rates(G, M, N):

    # equation 1
    μ1 = (V_G * G) / (K_G + G) # h^-1

    # equation 2
    μ2 = (V_M * M) / (K_M + M) * K_G_prime / (K_G_prime + G) # h^-1

    # equation 3
    μ3 = (V_N * N) / (K_N + N) * K_G_prime / (K_G_prime + G) * K_M_prime / (K_M_prime + M) # h^-1

    return μ1, μ2, μ3


# =============
# ODEs
# =============

def fermentation_odes(t, y):

    G, M, N, X, E, VDK, EA = y

    μ1, μ2, μ3 = specific_rates(G, M, N)

    # equation 4
    dG_dt = -μ1 * X # mol/(m^3*h)

    # equation 5
    dM_dt = -μ2 * X # mol/(m^3*h)

    # equation 6
    dN_dt = -μ3 * X # mol/(m^3*h)

    # equation 7
    growth_rate = (μ1 + μ2 + μ3) * X * (1 - X / X_max) # million cells/(mL*h)
    death_rate = k_d * X
    dX_dt = growth_rate - death_rate

    # equation 8
    dE_dt = Y_EG * (-dG_dt) + Y_EM * (-dM_dt) + Y_EN * (-dN_dt) # mol/(m^3*h)

    # equation 9
    dVDK_dt = k_f * X - k_r * VDK # ppm/h

    # equation 10
    dEA_dt = k_EA * X * E - k_hyd * EA # ppm/h

    return np.array([dG_dt, dM_dt, dN_dt, dX_dt, dE_dt, dVDK_dt, dEA_dt])


# =============
# RK4 Solver
# =============

def rk4_step(f, t, y, dt):

    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt * k1 / 2)
    k3 = f(t + dt / 2, y + dt * k2 / 2)
    k4 = f(t + dt, y + dt * k3)
    return y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


# =============
# Solve ODEs
# =============

# time array
t = np.linspace(t0, tf, n_points)
dt = t[1] - t[0]

# solution array
solution = np.zeros((n_points, 7))
solution[0] = [G0, M0, N0, X0, E0, VDK0, EA0]

for i in range(n_points - 1):
    solution[i + 1] = rk4_step(fermentation_odes, t[i], solution[i], dt)

glucose = solution[:, 0]
maltose = solution[:, 1]
maltotriose = solution[:, 2]
biomass = solution[:, 3]
ethanol = solution[:, 4]
vdks = solution[:, 5]
ethyl_acetate = solution[:, 6]

glucose_rate = np.zeros(n_points)
maltose_rate = np.zeros(n_points)
maltotriose_rate = np.zeros(n_points)

for i in range(n_points):
    glucose_rate[i], maltose_rate[i], maltotriose_rate[i] = specific_rates(
        glucose[i], maltose[i], maltotriose[i]
    )

# =============
# Save to CSV
# =============

with open('isothermal_fermentation_results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Time_h', 'Glucose', 'Maltose', 'Maltotriose', 'Biomass',
                     'Ethanol', 'VDKs', 'EthylAcetate', 'GlucoseRate',
                     'MaltoseRate', 'MaltotrioseRate'])

    for i in range(n_points):
        writer.writerow([
            t[i],
            glucose[i],
            maltose[i],
            maltotriose[i],
            biomass[i],
            ethanol[i],
            vdks[i],
            ethyl_acetate[i],
            glucose_rate[i],
            maltose_rate[i],
            maltotriose_rate[i]
        ])

# ===============
# Modified BFI needed for Sensitivity Analysis
# ===============

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
        VDK0,
        EA0,
    ])

    # time array
    t = np.linspace(t0, tf, n_points)
    dt = t[1] - t[0]

    solution = np.zeros((n_points, 7))
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

# ===============
# Baseline Operating Point
# ===============

G0_base = G0
M0_base = M0
N0_base = N0
T0_base = T0   # even if unused

baseline = run_fermentation(G0_base, M0_base, N0_base)

# ===============
# 3) Reusable Sensitivity Engine
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

# ===============
# Tornado plots (per output)
# ===============

def tornado_plot(sens_dict, title):
    labels = list(sens_dict.keys())
    values = [sens_dict[k] for k in labels]

    order = sorted(range(len(values)), key=lambda i: abs(values[i]))
    labels = [labels[i] for i in order]
    values = [values[i] for i in order]

    plt.figure()
    plt.barh(labels, values)
    plt.axvline(0, linewidth=0.8)
    plt.xlabel("Normalized sensitivity")
    plt.title(title)
    plt.tight_layout()
    plt.show()

for output, sens in sensitivities.items():
    tornado_plot(
        sens,
        f"Tornado plot – {output.replace('_', ' ')}"
    )

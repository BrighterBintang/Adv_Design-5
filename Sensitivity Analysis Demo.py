import numpy as np
import matplotlib.pyplot as plt

"""
SENSITIVITY ANALYSIS with Dummy Fermentation Model
(soon to be adapted to Kostas's ISOTHERMAL FERMENTATION MODEL on Thursday, 5 Feb 2026)
(based on Dimitrios's email saying that we may be going "too far", this sensitivity analysis may be scrapped)

Code Logic:
1) Dummy Fermentation Model
2) Baseline Operating Point
3) Reusable Local Sensitivity Engine
4) Computing the Sensitivity Analysis
5) One-At-a-Time Response Curves (Validating Linearity)
6) Tornado plots (per output)
"""

# ===============
# 1) Dummy Fermentation Model
# ===============

def fermentation_model(
    m_s0,
    T0,
    t_end=120,
    dt=0.1,
):
    """
    Simple dynamic beer fermentation model.

    Parameters
    ----------
    m_s0 : float
        Initial sugar mass (kg)
    T0 : float
        Initial temperature (°C)

    Returns
    -------
    dict with:
        time
        sugar
        ethanol
        max_rate
        fermentation_time
    """

    # --- constants ---
    Y_e = 0.51          # ethanol yield (kg/kg sugar)
    mu_max_20 = 0.4     # 1/h at 20°C
    Ks = 5.0            # Monod constant (kg)
    Ea = 0.07           # temperature sensitivity

    # temperature-adjusted growth rate
    mu_max = mu_max_20 * np.exp(Ea * (T0 - 20))

    # time grid
    t = np.arange(0, t_end, dt)

    sugar = np.zeros_like(t)
    ethanol = np.zeros_like(t)

    sugar[0] = m_s0

    rates = []

    for i in range(1, len(t)):
        mu = mu_max * sugar[i-1] / (Ks + sugar[i-1])
        dsdt = -mu * sugar[i-1]

        sugar[i] = max(sugar[i-1] + dsdt * dt, 0)
        ethanol[i] = Y_e * (m_s0 - sugar[i])

        rates.append(abs(dsdt))

        if sugar[i] < 0.01 * m_s0:
            sugar[i:] = sugar[i]
            ethanol[i:] = ethanol[i]
            break

    fermentation_time = t[i]
    max_rate = max(rates)

    return {
        "time": t,
        "sugar": sugar,
        "ethanol": ethanol,
        "final_ethanol": ethanol[i],
        "fermentation_time": fermentation_time,
        "max_rate": max_rate,
    }

# ===============
# 2) Baseline Operating Point
# ===============

m_s0 = 100.0   # Initial mass of sugars (kg)
T0 = 20.0      # Initial temperature (°C)

baseline = fermentation_model(m_s0, T0)

# ===============
# 3) Reusable Local Sensitivity Engine
# ===============

def local_sensitivity(
    model,
    m_s0,
    T0,
    dm_frac=0.05,
    dT=2.0,
):
    # baseline
    base = model(m_s0, T0)

    outputs = {
        "final_ethanol": base["final_ethanol"],
        "fermentation_time": base["fermentation_time"],
        "max_rate": base["max_rate"],
    }

    # sugar perturbation
    dm = dm_frac * m_s0
    plus_m = model(m_s0 + dm, T0)
    minus_m = model(m_s0 - dm, T0)

    # temperature perturbation
    plus_T = model(m_s0, T0 + dT)
    minus_T = model(m_s0, T0 - dT)

    results = {}

    for key, y0 in outputs.items():
        dy_dm = (plus_m[key] - minus_m[key]) / (2 * dm)
        dy_dT = (plus_T[key] - minus_T[key]) / (2 * dT)

        S_m = dy_dm * (m_s0 / y0)
        S_T = dy_dT * (T0 / y0)

        results[key] = {
            "S_sugar": S_m,
            "S_temp": S_T,
        }

    return results

# ===============
# 4) Computing the Sensitivity Analysis
# ===============

sens = local_sensitivity(fermentation_model, m_s0, T0)

for output, vals in sens.items():
    print(f"\n{output}")
    print(f"  Sugar sensitivity: {vals['S_sugar']:.3f}")
    print(f"  Temp sensitivity : {vals['S_temp']:.3f}")

# ===============
# 5) One-At-a-Time Response Curves (Validating Linearity)
# ===============

# Sugar sweep
m_vals = np.linspace(0.95*m_s0, 1.05*m_s0, 30)
ethanol_m = [
    fermentation_model(m, T0)["final_ethanol"] for m in m_vals
]

# Temperature sweep
T_vals = np.linspace(T0 - 2, T0 + 2, 30)
ethanol_T = [
    fermentation_model(m_s0, T)["final_ethanol"] for T in T_vals
]

plt.figure()
plt.plot(m_vals, ethanol_m)
plt.xlabel("Initial sugar mass (kg)")
plt.ylabel("Final ethanol (kg)")
plt.title("Ethanol response to sugar mass")
plt.show()

plt.figure()
plt.plot(T_vals, ethanol_T)
plt.xlabel("Initial temperature (°C)")
plt.ylabel("Final ethanol (kg)")
plt.title("Ethanol response to temperature")
plt.show()

# ===============
# 6) Tornado plots (per output)
# ===============

def tornado_plot(sens_dict, title):
    labels = ["Sugar mass", "Temperature"]
    values = [
        abs(sens_dict["S_sugar"]),
        abs(sens_dict["S_temp"]),
    ]

    plt.figure()
    plt.barh(labels, values)
    plt.xlabel("Absolute normalized sensitivity")
    plt.title(title)
    plt.show()

for output, vals in sens.items():
    tornado_plot(vals, f"Sensitivity of {output}")

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

x_data = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 
                   0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0])
y_data = np.array([0.0, 0.134, 0.23, 0.304, 0.365, 0.418, 0.496, 0.529, 0.573, 0.613, 
                   0.652, 0.697, 0.753, 0.818, 0.895, 0.946, 1.0])
get_y_theo = interp1d(x_data, y_data, kind='linear', fill_value="extrapolate")

M_eth, M_eau = 46.07, 18.015
H_vap_eth, H_vap_eau = 38560, 40660
rho_eth, rho_eau = 0.789, 1.000

def mole_to_abv(y):
    m_eth = y * M_eth
    m_eau = (1-y) * M_eau
    v_eth = m_eth / rho_eth
    v_eau = m_eau / rho_eau
    return 100 * v_eth / (v_eth + v_eau)

def simuler_passe(vol_init, abv_init, puissance):
    vol_eth = vol_init * abv_init
    vol_eau = vol_init * (1-abv_init)
    moles_eth = (vol_eth * 789) / M_eth
    moles_eau = (vol_eau * 1000) / M_eau
    m_eth_coll = 0
    m_eau_coll = 0
    t = 0
    dt = 5.0
    abv_hist, temps = [], []
    
    while moles_eth + moles_eau > 1.0 and (moles_eth/(moles_eth+moles_eau)) > 0.005:
        x = moles_eth/(moles_eth+moles_eau)
        y = float(get_y_theo(x))
        H = y*H_vap_eth + (1-y)*H_vap_eau
        dV = (puissance * dt) / H
        
        m_eth_coll += dV * y
        m_eau_coll += dV * (1-y)
        moles_eth -= dV * y
        moles_eau -= dV * (1-y)
        if moles_eth < 0: moles_eth = 0
        
        abv_hist.append(mole_to_abv(y))
        temps.append(t/3600)
        t += dt
        
    m_tot = m_eth_coll + m_eau_coll
    if m_tot > 0:
        y_final = m_eth_coll / m_tot
        vol_final = (m_eth_coll*M_eth/789) + (m_eau_coll*M_eau/1000)
        return vol_final, mole_to_abv(y_final)/100, temps, abv_hist
    return 0, 0, [], []

# CASCADE
v1, a1, t1, y1 = simuler_passe(50.0, 0.10, 2500) # Passe 1
v2, a2, t2, y2 = simuler_passe(v1, a1, 1500)     # Passe 2
v3, a3, t3, y3 = simuler_passe(v2, a2, 800)      # Passe 3

plt.figure(figsize=(10,6))
plt.plot(t1, y1, 'b:', label=f'Passe 1 (Moy: {a1*100:.1f}%)')
plt.plot(t2, y2, 'g--', label=f'Passe 2 (Moy: {a2*100:.1f}%)')
plt.plot(t3, y3, 'r-', linewidth=2, label=f'Passe 3 (Moy: {a3*100:.1f}%)')
plt.title("Gain de Pureté par Distillations Successives")
plt.ylabel("% ABV")
plt.xlabel("Temps (h)")
plt.legend()
plt.grid(True)
plt.show()
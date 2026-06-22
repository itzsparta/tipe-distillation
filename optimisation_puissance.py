import numpy as np
from scipy.interpolate import interp1d

# --- DONNÉES ---
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

def lancer_simulation(puissance_watts):
    # Setup 15L à 40% (Spirit Run)
    vol_eth = 15.0 * 0.40
    vol_eau = 15.0 * 0.60
    moles_eth = (vol_eth * rho_eth * 1000) / M_eth
    moles_eau = (vol_eau * rho_eau * 1000) / M_eau
    x_cuve = moles_eth / (moles_eth + moles_eau)
    L_moles = moles_eth + moles_eau
    
    dt = 5.0
    t = 0.0
    
    # Modèle Efficacité (Perte si > 800W)
    efficacite = 1.0 - (puissance_watts - 800)*0.00015 if puissance_watts > 800 else 1.0
    if efficacite < 0.6: efficacite = 0.6
    
    y_reel_list, vol_cumule = [], []
    vol_prod = 0.0
    
    while L_moles > 1.0 and x_cuve > 0.001:
        y_theo = float(get_y_theo(x_cuve))
        y_reel = y_theo * efficacite
        H_vap = y_reel * H_vap_eth + (1 - y_reel) * H_vap_eau
        dV = (puissance_watts * dt) / H_vap
        
        M_dist = y_reel * M_eth + (1-y_reel) * M_eau
        rho_dist = y_reel * 0.79 + (1-y_reel) * 1.0
        vol_prod += (dV * M_dist / 1000) / rho_dist
        
        moles_eth -= dV * y_reel
        moles_eau -= dV * (1-y_reel)
        if moles_eth < 0: moles_eth = 0
        L_moles = moles_eth + moles_eau
        x_cuve = moles_eth / L_moles
        
        y_reel_list.append(y_reel)
        vol_cumule.append(vol_prod)
        t += dt
        if t > 10*3600: break

    # Analyse Cœur
    abv_hist = np.array([mole_to_abv(y) for y in y_reel_list])
    idx_tetes = np.searchsorted(vol_cumule, 0.200)
    idx_queues = np.where(abv_hist[idx_tetes:] < 55.0)[0]
    idx_fin = idx_tetes + idx_queues[0] if len(idx_queues)>0 else len(abv_hist)
    
    if idx_fin > idx_tetes:
        return {
            "P": puissance_watts, "Eff": efficacite*100,
            "ABV": np.mean(abv_hist[idx_tetes:idx_fin]),
            "Vol": vol_cumule[idx_fin] - vol_cumule[idx_tetes],
            "Temps": t/3600,
            "E": (puissance_watts * (t*(idx_fin-idx_tetes)/len(abv_hist)))/3.6e6
        }
    return None

# Affichage Tableau
print(f"{'Puissance':^10} | {'Efficacité':^10} | {'Pureté Cœur':^12} | {'Vol. Cœur':^10} | {'Durée':^10}")
print("-" * 65)
for P in [800, 1200, 1600, 2000, 2400]:
    r = lancer_simulation(P)
    print(f"{r['P']:^10} | {r['Eff']:^10.1f} | {r['ABV']:^12.2f} | {r['Vol']:^10.2f} | {r['Temps']:^10.2f}")
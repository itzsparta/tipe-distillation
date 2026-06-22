import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- 1. DONNÉES & INTERPOLATION ---
x_data = np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.40, 
                   0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0])
y_data = np.array([0.0, 0.134, 0.23, 0.304, 0.365, 0.418, 0.496, 0.529, 0.573, 0.613, 
                   0.652, 0.697, 0.753, 0.818, 0.895, 0.946, 1.0])
t_data = np.array([100.0, 96.4, 93.5, 91.2, 89.8, 88.3, 86.0, 84.1, 81.9, 80.7, 
                   79.8, 79.2, 78.7, 78.3, 78.15, 78.15, 78.3]) 

get_y_vapeur = interp1d(x_data, y_data, kind='linear', fill_value="extrapolate")
get_T_ebullition = interp1d(x_data, t_data, kind='linear', fill_value="extrapolate")

M_eth, M_eau = 46.07, 18.015
H_vap_eth, H_vap_eau = 38560, 40660
rho_eth, rho_eau = 0.789, 1.000

# --- 2. CONFIGURATION (15L Moût Neutre) ---
volume_mout_L = 15.0  
degre_alcool_vol = 0.12 # 12%

vol_eth = volume_mout_L * degre_alcool_vol
vol_eau = volume_mout_L * (1.0 - degre_alcool_vol)
moles_eth = (vol_eth * rho_eth * 1000) / M_eth
moles_eau = (vol_eau * rho_eau * 1000) / M_eau
L_total_moles = moles_eth + moles_eau
x_eth_cuve = moles_eth / L_total_moles

# --- 3. PARAMÈTRES CHAUFFE DYNAMIQUE ---
Q_base = 800.0   
dt = 5.0
temps_max = 5 * 3600 
increment_watts = 100.0
intervalle_increment_s = 20 * 60 

# --- 4. SIMULATION ---
temps, y_sortie, debit_sortie = [], [], []
volume_cumule_distillat_L, Q_history = [], []
vol_total_produit = 0.0
t = 0.0

while t < temps_max and L_total_moles > 1.0 and x_eth_cuve > 0.001:
    # Chauffe escalier
    nb_paliers = int(t / intervalle_increment_s)
    Q_actuel = Q_base + (nb_paliers * increment_watts)
    if Q_actuel > 2000: Q_actuel = 2000 

    try:
        y_eth_vapeur = float(get_y_vapeur(x_eth_cuve))
    except: break

    H_vap_mix = y_eth_vapeur * H_vap_eth + (1 - y_eth_vapeur) * H_vap_eau
    dV_moles = (Q_actuel * dt) / H_vap_mix
    
    # Volume produit
    M_dist = y_eth_vapeur * M_eth + (1-y_eth_vapeur) * M_eau
    rho_dist = y_eth_vapeur * 0.79 + (1-y_eth_vapeur) * 1.0
    vol_produit_L = (dV_moles * M_dist / 1000) / rho_dist
    vol_total_produit += vol_produit_L

    # Mise à jour Cuve
    moles_eth -= dV_moles * y_eth_vapeur
    moles_eau -= dV_moles * (1.0 - y_eth_vapeur)
    if moles_eth < 0: moles_eth = 0
    L_total_moles = moles_eth + moles_eau
    x_eth_cuve = moles_eth / L_total_moles
    
    temps.append(t / 3600.0)
    y_sortie.append(y_eth_vapeur)
    debit_sortie.append(dV_moles / dt)
    volume_cumule_distillat_L.append(vol_total_produit)
    Q_history.append(Q_actuel)
    t += dt

# --- 5. GRAPHIQUES ---
def mole_to_abv(y):
    m_eth = y * M_eth
    m_eau = (1-y) * M_eau
    v_eth = m_eth / rho_eth
    v_eau = m_eau / rho_eau
    return 100 * v_eth / (v_eth + v_eau)

distillat_abv = np.array([mole_to_abv(y) for y in y_sortie])
vol_cumule = np.array(volume_cumule_distillat_L)
temps_arr = np.array(temps)

# Indices Coupes (Têtes 200ml / Queues < 40%)
idx_fin_tetes = np.searchsorted(vol_cumule, 0.2) 
idx_queues = np.where(distillat_abv[idx_fin_tetes:] < 40.0)[0]
idx_debut_queues = idx_fin_tetes + idx_queues[0] if len(idx_queues) > 0 else len(distillat_abv)

plt.figure(figsize=(12, 10))
plt.subplot(2, 1, 1)
plt.plot(temps_arr[:idx_fin_tetes], distillat_abv[:idx_fin_tetes], 'r', label='Têtes', linewidth=2)
plt.plot(temps_arr[idx_fin_tetes:idx_debut_queues], distillat_abv[idx_fin_tetes:idx_debut_queues], 'g', label='Cœur', linewidth=2)
plt.plot(temps_arr[idx_debut_queues:], distillat_abv[idx_debut_queues:], 'gray', label='Queues', linestyle='--')
plt.ylabel('% ABV')
plt.title("Simulation Coupes (Têtes/Cœur/Queues)")
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
ax1 = plt.gca()
ax1.plot(temps, debit_sortie, 'k-', label='Débit')
ax1.set_ylabel('Débit (mol/s)')
ax2 = ax1.twinx()
ax2.plot(temps, Q_history, 'm:', label='Puissance')
ax2.set_ylabel('Puissance (W)', color='m')
plt.title("Débit et Puissance")
plt.grid(True)
plt.show()
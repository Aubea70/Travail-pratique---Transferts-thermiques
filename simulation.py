import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as date
import pandas as pd


#lire les données
data = pd.read_csv("Dataset of weighing station temperature measurements.csv",delimiter=";")
time_str = data["Time"]
temp = data["Outdoor temperature [deg. C]"]
humid = data["Outdoor relative humidity [%]"]
low = data.loc[:,"T[degC]-Low-S1":"T[degC]-Low-S29"].to_numpy()
mid = data.loc[:,"T[degC]-Mid-S1":"T[degC]-Mid-S29"].to_numpy()
top = data.loc[:,"T[degC]-Top-S1":"T[degC]-Top-S29"].to_numpy()

format = "%Y-%m-%d %H:%M"
time = np.array([date.strptime(t, format) for t in time_str])


# -----------------------
# PARAMÈTRES DU MODÈLE
# -----------------------

dt = 120                 # pas de temps (s) = 2 min
N_jours = 120            # durée de la simu (jours) ~ 4 mois
N = int(N_jours * 24 * 3600 / dt)

# Propriétés des matériaux et dimensions du puit

rho_air, cp_air = 1.2, 1000        
rho_beton, cp_beton, k_beton = 2350, 1000, 1.4   

L, W, H = 26.1, 3.7, 1.7

# Air
V_air = L * W * H   # Vérifier si on enlève les colonnes de béton du volume
C_air = rho_air * cp_air * V_air
h_int = 10.45       # Mettre les sources des constantes dans le rapport
h_ext = 27.8
A_craques = 2e-3 * (W*7 + L*2)
vit_air = 12    # vitesse moyenne de l'air en hiver[m/s]
débit = rho_air*vit_air*A_craques/2     # divisé par 2 car moitié des craques entrée et l'autre moitié sortie

# Dalle plafond
e_dalle = 0.4       # On assume l'épaisseur plafond = épaisseur béton. Vérifier si vrai.
A_dalle = L * W     # Hypothèse : craques entre les plaques négligeables
V_plaf = A_dalle * e_dalle
C_plafond = rho_beton * cp_beton * V_plaf

# Dalle plancher
C_plancher = C_plafond

# Mur
A_mur = 2*(L + W) * H
e_mur = 0.4
V_mur = A_mur * e_mur
C_mur = rho_beton * cp_beton * V_mur

# Isolant
e_iso = 0.1
k_iso = 0.035


# Résistance conduction : On néglige les arrêtes et les coins (ou si on prend un facteur de forme?)
R_dalle = e_dalle / (k_beton * A_dalle)    
R_mur_cond = e_mur / (k_beton * A_mur)
R_iso_dalle = e_iso / (k_iso * A_dalle)



# Convection
R_air_plaf_conv = 1 / (h_int * A_dalle) # Même résistance pour le plancher
R_plaf_ext_conv = 1 / (h_ext * A_dalle)
R_air_mur_conv  = 1 / (h_int * A_mur)


# Plafond
R_air_plaf = R_air_plaf_conv
R_plaf_ext = R_plaf_ext_conv + R_dalle

# Mur

R_air_mur = R_air_mur_conv + R_mur_cond
R_iso_mur= e_iso / (k_iso * A_mur) + R_mur_cond

# Plancher

R_air_plan = R_air_plaf_conv
R_plan_iso = R_iso_dalle + R_dalle



# --------------------------
# CALCUL DES FLUX DE CHALEUR
# --------------------------

matrice  = [    #pour T_in, T_plaf, T_mur, T_sol, T_ex, T_surfterre
            [(C_air/dt + 2/R_air_plaf_conv + 1/R_air_mur_conv + débit*cp_air), -1/R_air_plaf_conv, -1/R_air_mur_conv, -1/R_air_plaf_conv, -débit*cp_air, 0],
            [-1/R_air_plaf_conv, (C_plafond/dt + 1/R_air_plaf_conv + 1/R_plaf_ext), 0, 0, -1/R_plaf_ext, 0],
            [-1/R_air_mur_conv, 0, (C_mur/dt + 1/R_air_mur_conv + 1/R_iso_mur), 0, 0, -1/(2*R_iso_mur)],
            [-1/R_air_plaf_conv, 0, 0, (C_plancher/dt + 1/R_air_plaf_conv + 1/R_plan_iso), 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, -1/R_iso_mur, 0, 0, (1/(2*R_iso_mur) + 1/R_plaf_ext_conv)],
            ]

# M*T(t) = T(t-dt) -> T(t) = M**(-1)*T(t-dt)
mat_inv = np.linalg.inv(matrice)

T_in = np.zeros(len(time))
T_plaf = np.zeros(len(time))
T_mur = np.zeros(len(time))
T_sol = np.zeros(len(time))
T_ex = temp.to_numpy()
T_surfterre = np.zeros(len(time))
T_terre = 8     # Température constante à 8 deg C sous le puits
q_aero = np.zeros(len(time))

# Conditions initiales:
T_in[0] = np.nanmean(data.loc[0,"T[degC]-Low-S1":"T[degC]-Top-S29"].to_numpy())      # Température air interne initiale (moyenne dans le puits initiale)
T_plaf[0] = np.nanmean(top[0,:])
T_mur[0] = np.nanmean(mid[0,:])
T_sol[0] = np.nanmean(low[0,:])
T_surfterre[0] = T_ex[0]     # Assume surface de la terre température proche de l'air externe
if T_ex[0] < 0:
    q_aero[0] = 60e3

# Calcul
for i in range(1,len(time)):
    if T_ex[i] < 0:
        q_aero[i] = 60e3    #puissance*dt

    output = np.matmul(mat_inv, [C_air * T_in[i-1]/dt + q_aero[i-1],
                                 C_plafond * T_plaf[i-1]/dt,
                                 C_mur * T_mur[i-1]/dt + 4/R_iso_mur,
                                 C_plancher * T_sol[i-1]/dt + 8/R_plan_iso,
                                 T_ex[i-1],
                                 -4/R_iso_mur])
    T_in[i], T_plaf[i], T_mur[i], T_sol[i], _, T_surfterre[i] = output


sim_data = pd.DataFrame({"Time":time, "T_in":T_in, "T_plaf":T_plaf, "T_mur":T_mur, "T_sol":T_sol, "T_surfterre":T_surfterre})
sim_data.to_csv("Simulated_temperatures.csv", index=False)

énergie = sum(q_aero*dt)

print("Énergie totale consommée :", (énergie/1000)/360, "kWh")

plt.plot(time, T_in, label="Simulée")
plt.plot(time, np.nanmean(data.loc[:,"T[degC]-Low-S1":"T[degC]-Top-S29"].to_numpy(), axis=1), label="Réelle")
plt.xlabel("Temps")
plt.ylabel("Température [°C]")
plt.legend()
plt.show()
import numpy as np
import matplotlib.pyplot as plt

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
débit = ...

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



# Conditions initiales:
data = ...      # Températures externes à simuler
T_ex = data[0]  # Température de l'air externe
T_gc = 8        # Température constante à 8 deg C sous le puits
T_gs = T_ex     # Température de la surface de la terre exposée à l'air extérieur
T_in = ...      # Température air interne initiale (prendre 1ère valeur des données ou imposer la température seuil de 3 deg?)


# --------------------------
# CALCUL DES FLUX DE CHALEUR
# --------------------------

matrice  = [    #pour T_a, T_p, T_m, T_s, T_inf, T_gs, 1
            [(C_air/dt + 2/R_air_plaf_conv + 1/R_air_mur_conv + débit*cp_air), -1/R_air_plaf_conv, -1/R_air_mur_conv, -1/R_air_plaf_conv, 0, 0, 0],
            [-1/R_air_plaf_conv, (C_plafond/dt + 1/R_air_plaf_conv + 1/R_plaf_ext), 0, 0, -1/R_plaf_ext, 0, 0],
            [-1/R_air_mur_conv, 0, (C_mur/dt + 1/R_air_mur_conv + 1/R_iso_mur), 0, 0, -1/(2*R_iso_mur), -4/R_iso_mur],
            [-1/R_air_plaf_conv, 0, 0, (C_plancher/dt + 1/R_air_plaf_conv + 1/R_plan_iso), 0, 0, -8/R_plan_iso],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 0, -1/R_iso_mur, 0, 0, (1/(2*R_iso_mur) + 1/R_plaf_ext_conv), 4/R_iso_mur],
            [0, 0, 0, 0, 0, 0, 1]
            ]

# M*T(t) = T(t-dt) -> T(t) = M**(-1)*T(t-dt)
mat_inv = np.linalg.inv(matrice)
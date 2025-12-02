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
V_air = L * W * H
C_air = rho_air * cp_air * V_air
h_int = 10.45
h_ext = 27.8

# Dalle plafond
e_dalle = 0.4
A_dalle = L * W
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
k_iso = 0.0375


# Résistance conduction
R_dalle = e_dalle / (k_beton * A_dalle)    
R_mur_cond = e_mur / (k_beton * A_mur)
R_iso_dalle = e_iso / (k_iso * A_dalle)


# Convection
R_air_plaf_conv = 1 / (h_int * A_dalle) # Même résistance pour le plancher
R_plaf_ext_conv = 1 / (h_ext * A_dalle)
R_air_mur_conv  = 1 / (h_int * A_mur)


# Plafond
R_air_plaf = R_dalle + R_air_plaf_conv
R_plaf_ext = R_plaf_ext_conv

# Mur

R_air_mur = R_air_mur_conv + R_mur_cond
R_iso_mur= e_iso / (k_iso * A_mur)

# Plancher

R_air_plan = R_air_plaf_conv + R_dalle
R_plan_iso = R_iso_mur

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime as date
import numpy as np

# --- Load CSV ---
df = pd.read_csv(r"C:\Users\laura\Documents\Education\01 - Université\Génie Physique\A25\TT\tp\Dataset of weighing station temperature measurements.csv", sep=";")

# --- Columns to average ---
average_low_P1 = [f"T[degC]-Low-S{i}" for i in range(1,6)]
average_mid_P1 = [f"T[degC]-Mid-S{i}" for i in range(1,6)]
average_top_P1 = [f"T[degC]-Top-S{i}" for i in range(1,6)]

average_low_P2 = [f"T[degC]-Low-S{i}" for i in range(5,10)]
average_mid_P2 = [f"T[degC]-Mid-S{i}" for i in range(5,10)]
average_top_P2 = [f"T[degC]-Top-S{i}" for i in range(5,10)]

average_low_P3 = [f"T[degC]-Low-S{i}" for i in range(9,14)]
average_mid_P3 = [f"T[degC]-Mid-S{i}" for i in range(9,14)]
average_top_P3 = [f"T[degC]-Top-S{i}" for i in range(9,14)]

average_low_P4 = [f"T[degC]-Low-S{i}" for i in range(13,20)]
average_mid_P4 = [f"T[degC]-Mid-S{i}" for i in range(13,20)]
average_top_P4 = [f"T[degC]-Top-S{i}" for i in range(13,20)]

average_low_P5 = [f"T[degC]-Low-S{i}" for i in range(19,24)]
average_mid_P5 = [f"T[degC]-Mid-S{i}" for i in range(19,24)]
average_top_P5 = [f"T[degC]-Top-S{i}" for i in range(19,24)]

average_low_P6 = [f"T[degC]-Low-S{i}" for i in range(23,30)]
average_mid_P6 = [f"T[degC]-Mid-S{i}" for i in range(23,30)]
average_top_P6 = [f"T[degC]-Top-S{i}" for i in range(23,30)]

# --- Compute row-wise mean ---
df["row_mean_low_P1"] = df[average_low_P1].mean(axis=1)
df["row_mean_mid_P1"] = df[average_mid_P1].mean(axis=1)
df["row_mean_top_P1"] = df[average_top_P1].mean(axis=1)

df["row_mean_low_P2"] = df[average_low_P2].mean(axis=1)
df["row_mean_mid_P2"] = df[average_mid_P2].mean(axis=1)
df["row_mean_top_P2"] = df[average_top_P2].mean(axis=1)

df["row_mean_low_P3"] = df[average_low_P3].mean(axis=1)
df["row_mean_mid_P3"] = df[average_mid_P3].mean(axis=1)
df["row_mean_top_P3"] = df[average_top_P3].mean(axis=1)

df["row_mean_low_P4"] = df[average_low_P4].mean(axis=1)
df["row_mean_mid_P4"] = df[average_mid_P4].mean(axis=1)
df["row_mean_top_P4"] = df[average_top_P4].mean(axis=1)

df["row_mean_low_P5"] = df[average_low_P5].mean(axis=1)
df["row_mean_mid_P5"] = df[average_mid_P5].mean(axis=1)
df["row_mean_top_P5"] = df[average_top_P5].mean(axis=1)

df["row_mean_low_P6"] = df[average_low_P6].mean(axis=1)
df["row_mean_mid_P6"] = df[average_mid_P6].mean(axis=1)
df["row_mean_top_P6"] = df[average_top_P6].mean(axis=1)

# --- Compute cell mean ---
df["P1_mean"] = (df["row_mean_low_P1"] + df["row_mean_mid_P1"] + df["row_mean_top_P1"])/3

df["P2_mean"] = (df["row_mean_low_P2"] + df["row_mean_mid_P2"] + df["row_mean_top_P2"])/3

df["P3_mean"] = (df["row_mean_low_P3"] + df["row_mean_mid_P3"] + df["row_mean_top_P3"])/3

df["P4_mean"] = (df["row_mean_low_P4"] + df["row_mean_mid_P4"] + df["row_mean_top_P4"])/3

df["P5_mean"] = (df["row_mean_low_P5"] + df["row_mean_mid_P5"] + df["row_mean_top_P5"])/3

df["P6_mean"] = (df["row_mean_low_P6"] + df["row_mean_mid_P6"] + df["row_mean_top_P6"])/3

# --- Convert Time to datetime ---
df['Time'] = pd.to_datetime(df['Time'], errors='coerce')

out_temp = "Outdoor temperature [deg. C]"


# --- Downsample for plotting ---
df_sample = df.iloc[::50].copy()


# Create 2 vertically stacked subplots that share the x-axis
fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(12, 8), sharex=True,
    gridspec_kw={'height_ratios': [3, 1, 1]}  # top graph bigger
)

# ----------------------
# TOP GRAPH – TEMPERATURE
# ----------------------
data = pd.read_csv("Dataset of weighing station temperature measurements.csv",delimiter=";")
temp = data["Outdoor temperature [deg. C]"]
time_str = data["Time"]
format = "%Y-%m-%d %H:%M"
time = np.array([date.strptime(t, format) for t in time_str])
with np.errstate(all='ignore'):
    mean_values = np.nanmean(data.loc[:, "T[degC]-Low-S1":"T[degC]-Top-S29"].to_numpy(), axis=1)

data2 = pd.read_csv("intemperies_2min_complet.csv", delimiter=";")  # ou ',' selon le fichier
print(data2.columns)

time2 = pd.to_datetime(data2.iloc[:, 0])  # première colonne = temps
valeurs = data2.iloc[:, 1]

# Graphique original

ax1.plot(time, mean_values, label="Température moyenne intérieure du puits")
ax1.plot(time, temp, label="Température extérieure")
ax1.axhline(y=0, color="black", linestyle="--")
ax1.set_ylabel("Température [°C]")
ax1.legend()

# Nouveau graphique
ax2.plot(time2, valeurs, color="orange", label="Précipitations")
ax2.set_xlabel("Temps")
ax2.set_ylabel("Quantité [mm]")
ax2.legend()

ax3.plot(df_sample["Time"], df_sample["Outdoor relative humidity [%]"], color='purple', label="Taux d'humidité")

ax3.set_ylabel("Humidité relative extérieure [%]")
ax3.set_xlabel("Temps")
ax3.grid(True)
ax3.legend()


# ----------------------------
# Auto-format x-axis & layout
# ----------------------------
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})
fig.autofmt_xdate()
plt.tight_layout()

# --- Save figure ---
plt.savefig(r"C:\Users\laura\Documents\Education\01 - Université\Génie Physique\A25\TT\tp\cell_plot.png", dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# --- Load CSV ---
df = pd.read_csv(r"C:\Users\laura\Documents\Education\01 - Université\Génie Physique\A25\TT\tp\Dataset of weighing station temperature measurements.csv", sep=";")

# --- Columns to average ---
columns_to_average_low = [f"T[degC]-Low-S{i}" for i in range(1,30)]
columns_to_average_mid = [f"T[degC]-Mid-S{i}" for i in range(1,30)]
columns_to_average_top = [f"T[degC]-Top-S{i}" for i in range(1,30)]

# --- Compute row-wise mean ---
df["row_mean_low"] = df[columns_to_average_low].mean(axis=1)
df["row_mean_mid"] = df[columns_to_average_mid].mean(axis=1)
df["row_mean_top"] = df[columns_to_average_top].mean(axis=1)

# --- Compute differences ---
df["diff_top_mid"] = df["row_mean_top"] - df["row_mean_mid"]
df["diff_mid_low"] = df["row_mean_mid"] - df["row_mean_low"]

# --- Compute overall averages ---
avg_diff_top_mid = df["diff_top_mid"].mean()
avg_diff_mid_low = df["diff_mid_low"].mean()

# --- Print results ---
print(f"Average difference Top-Mid: {avg_diff_top_mid:.2f} °C")
print(f"Average difference Mid-Low: {avg_diff_mid_low:.2f} °C")

# --- Convert Time to datetime ---
df['Time'] = pd.to_datetime(df['Time'], errors='coerce')

# --- Downsample for plotting ---
df_sample = df.iloc[::50].copy()

# --- Plot all three means on the same graph ---
plt.figure(figsize=(12,6))
plt.plot(df_sample['Time'], df_sample["row_mean_low"], color='blue', marker='o', markersize=3, linestyle='-', label='Capteurs bas')
plt.plot(df_sample['Time'], df_sample["row_mean_mid"], color='green', marker='s', markersize=3, linestyle='-', label='Capteurs mi-hauteur')
plt.plot(df_sample['Time'], df_sample["row_mean_top"], color='red', marker='^', markersize=3, linestyle='-', label='Capteurs hauts')

plt.xlabel("Temps", fontsize=14)
plt.ylabel("Température (°C)", fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.legend(fontsize=12)
plt.grid(True)
plt.gcf().autofmt_xdate()  # format dates nicely
plt.tight_layout()

# --- Save figure ---
plt.savefig(r"C:\Users\laura\Documents\Education\01 - Université\Génie Physique\A25\TT\tp\line_plot.png", dpi=300)
plt.show()

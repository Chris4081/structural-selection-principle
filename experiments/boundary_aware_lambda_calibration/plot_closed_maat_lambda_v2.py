import json
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

OUTDIR = "plots"
os.makedirs(OUTDIR, exist_ok=True)

# --- Load data ---
df = pd.read_csv("maat_defects_fused.csv")

# --- Closed MAAT lambda result ---
with open("closed_maat_lambda_fit_results.json", "r") as f:
    fit_result = json.load(f)

lambdas = fit_result["lambdas"]

labels = list(lambdas.keys())
values = np.array(list(lambdas.values()))
shares = values / values.sum()

# --- Figure 1: Lambda distribution ---
plt.figure(figsize=(7, 4.5))
plt.bar(labels, values)
plt.ylabel(r"Calibrated $\lambda_a$")
plt.title("Closed MAAT λ v2 — Fused Calibration")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig1_lambda_distribution.png"), dpi=300)
plt.close()

# --- Figure 2: Normalized shares ---
plt.figure(figsize=(7, 4.5))
plt.bar(labels, shares)
plt.ylabel("Normalized contribution")
plt.title("Normalized MAAT Weight Shares")
plt.ylim(0, max(shares) * 1.25)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig2_lambda_shares.png"), dpi=300)
plt.close()

# --- Figure 3: SAT vs CORE defect comparison ---
defect_cols = ["d_H", "d_B", "d_S", "d_V", "d_R"]

means = df.groupby("source")[defect_cols].mean().T
means.index = ["H", "B", "S", "V", "R"]

plt.figure(figsize=(8, 4.8))
x = np.arange(len(means.index))
width = 0.35

sources = list(means.columns)

for i, source in enumerate(sources):
    plt.bar(x + (i - 0.5) * width, means[source], width, label=source)

plt.xticks(x, means.index)
plt.ylabel("Mean defect")
plt.title("Mean Structural Defects: SAT vs MAAT-Core")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig3_sat_vs_core_defects.png"), dpi=300)
plt.close()

# --- Figure 4: structural score distribution ---
lambda_vec = np.array([
    lambdas["H"],
    lambdas["B"],
    lambdas["S"],
    lambdas["V"],
    lambdas["R"],
])

D = df[defect_cols].to_numpy()
df["F_lambda"] = D @ lambda_vec

plt.figure(figsize=(7, 4.5))
for source in df["source"].unique():
    subset = df[df["source"] == source]
    plt.hist(subset["F_lambda"], bins=40, alpha=0.55, label=source)

plt.xlabel(r"$F_\lambda = \sum_a \lambda_a d_a$")
plt.ylabel("Count")
plt.title("Structural Selection Energy Distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig4_structural_energy_distribution.png"), dpi=300)
plt.close()

# --- Figure 5: C_hat if available in separate SAT file ---
try:
    sat_raw = pd.read_csv("maat_defects_sat.csv")
    if "C_hat" in sat_raw.columns and "runtime" in sat_raw.columns:
        Dsat = sat_raw[defect_cols].to_numpy()
        sat_raw["F_lambda"] = Dsat @ lambda_vec

        plt.figure(figsize=(7, 4.5))
        plt.scatter(sat_raw["C_hat"], sat_raw["F_lambda"], s=10, alpha=0.45)
        plt.xlabel(r"$\hat{C}$")
        plt.ylabel(r"$F_\lambda$")
        plt.title(r"SAT: $\hat{C}$ vs Calibrated Structural Energy")
        plt.tight_layout()
        plt.savefig(os.path.join(OUTDIR, "fig5_c_hat_vs_flambda_sat.png"), dpi=300)
        plt.close()
except FileNotFoundError:
    pass

print("Figures written:")
print(os.path.join(OUTDIR, "fig1_lambda_distribution.png"))
print(os.path.join(OUTDIR, "fig2_lambda_shares.png"))
print(os.path.join(OUTDIR, "fig3_sat_vs_core_defects.png"))
print(os.path.join(OUTDIR, "fig4_structural_energy_distribution.png"))
print(os.path.join(OUTDIR, "fig5_c_hat_vs_flambda_sat.png"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# INPUT FILES
# ============================================================

FILE_1D = "cci_entropy_information_test.csv"
FILE_2D = "2d_cci_entropy_information_test.csv"

# Effective exponents from the papers
ALPHA_1D = 2.5
ALPHA_2D = 3.0   # lower bound within tested range

EPS = 1e-8

# ============================================================
# LOAD DATA
# ============================================================

df1 = pd.read_csv(FILE_1D)
df2 = pd.read_csv(FILE_2D)

# Required columns:
# mean_cci, mean_mi, mean_dS_pos
required_cols = ["mean_cci", "mean_mi", "mean_dS_pos"]

for col in required_cols:
    if col not in df1.columns:
        raise ValueError(f"1D file missing column: {col}")
    if col not in df2.columns:
        raise ValueError(f"2D file missing column: {col}")

# ============================================================
# BUILD SCALING RATIOS
# ============================================================

df1["ratio_best"] = df1["mean_dS_pos"] / ((df1["mean_mi"] + EPS) ** ALPHA_1D)
df2["ratio_best"] = df2["mean_dS_pos"] / ((df2["mean_mi"] + EPS) ** ALPHA_2D)

# Optional: log-x values for a cleaner comparison
df1["log_ratio_best"] = np.log10(df1["ratio_best"] + EPS)
df2["log_ratio_best"] = np.log10(df2["ratio_best"] + EPS)

# ============================================================
# FIGURE 2: SIDE-BY-SIDE SCALING COMPARISON
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)

# ------------------------------------------------------------
# Panel (a): 1D
# ------------------------------------------------------------
ax = axes[0]
ax.scatter(df1["ratio_best"], df1["mean_cci"], alpha=0.8)
ax.set_xlabel(r"$R_{\alpha_{1D}}=\overline{\dot S}_{cg}^{(+)}/(\overline{I}_{nn}+\varepsilon)^{2.5}$")
ax.set_ylabel("mean CCI")
ax.set_title(r"(a) 1D scaling relation")
ax.grid(alpha=0.3)

ax.text(
    0.05, 0.95,
    "1D\n" + r"$\alpha \approx 2.5$" + "\n" + r"$r_s \approx 0.734$",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

# ------------------------------------------------------------
# Panel (b): 2D
# ------------------------------------------------------------
ax = axes[1]
ax.scatter(df2["ratio_best"], df2["mean_cci"], alpha=0.8)
ax.set_xlabel(r"$R_{\alpha_{2D}}=\overline{\dot S}_{cg}^{(+)}/(\overline{I}_{nn}+\varepsilon)^{3.0}$")
ax.set_title(r"(b) 2D scaling relation")
ax.grid(alpha=0.3)

ax.text(
    0.05, 0.95,
    "2D\n" + r"$\alpha \gtrsim 3.0$" + "\n" + r"$r_s \approx 0.870$" + "\nwithin tested range",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.tight_layout()
plt.savefig("figure2_scaling_comparison.png", dpi=220)
plt.show()

# ============================================================
# OPTIONAL VERSION: LOG-X AXIS
# Uncomment if the x-ranges are too different
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)

# 1D log version
ax = axes[0]
ax.scatter(df1["log_ratio_best"], df1["mean_cci"], alpha=0.8)
ax.set_xlabel(r"$\log_{10} R_{\alpha_{1D}}$")
ax.set_ylabel("mean CCI")
ax.set_title(r"(a) 1D scaling relation (log-x)")
ax.grid(alpha=0.3)

ax.text(
    0.05, 0.95,
    "1D\n" + r"$\alpha \approx 2.5$" + "\n" + r"$r_s \approx 0.734$",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

# 2D log version
ax = axes[1]
ax.scatter(df2["log_ratio_best"], df2["mean_cci"], alpha=0.8)
ax.set_xlabel(r"$\log_{10} R_{\alpha_{2D}}$")
ax.set_title(r"(b) 2D scaling relation (log-x)")
ax.grid(alpha=0.3)

ax.text(
    0.05, 0.95,
    "2D\n" + r"$\alpha \gtrsim 3.0$" + "\n" + r"$r_s \approx 0.870$" + "\nwithin tested range",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.tight_layout()
plt.savefig("figure2_scaling_comparison_logx.png", dpi=220)
plt.show()
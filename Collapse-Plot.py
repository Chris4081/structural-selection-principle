import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# INPUT FILES
# ============================================================

FILE_1D = "cci_entropy_information_test.csv"
FILE_2D = "2d_cci_entropy_information_test.csv"
FILE_3D = "3d_cci_entropy_information_test.csv"

EPS = 1e-8

# Effective exponents from the papers / current results
ALPHA_1D = 2.5
ALPHA_2D = 3.0      # lower bound within tested range
ALPHA_3D = 2.8      # representative point from plateau region

# ============================================================
# LOAD DATA
# ============================================================

df1 = pd.read_csv(FILE_1D)
df2 = pd.read_csv(FILE_2D)
df3 = pd.read_csv(FILE_3D)

required_cols = ["mean_cci", "mean_mi", "mean_dS_pos"]

for col in required_cols:
    if col not in df1.columns:
        raise ValueError(f"1D file missing column: {col}")
    if col not in df2.columns:
        raise ValueError(f"2D file missing column: {col}")
    if col not in df3.columns:
        raise ValueError(f"3D file missing column: {col}")

# ============================================================
# BUILD SCALING RATIOS
# ============================================================

df1["ratio"] = df1["mean_dS_pos"] / ((df1["mean_mi"] + EPS) ** ALPHA_1D)
df2["ratio"] = df2["mean_dS_pos"] / ((df2["mean_mi"] + EPS) ** ALPHA_2D)
df3["ratio"] = df3["mean_dS_pos"] / ((df3["mean_mi"] + EPS) ** ALPHA_3D)

df1["log_ratio"] = np.log10(df1["ratio"] + EPS)
df2["log_ratio"] = np.log10(df2["ratio"] + EPS)
df3["log_ratio"] = np.log10(df3["ratio"] + EPS)

df1["log_cci"] = np.log10(df1["mean_cci"] + EPS)
df2["log_cci"] = np.log10(df2["mean_cci"] + EPS)
df3["log_cci"] = np.log10(df3["mean_cci"] + EPS)

# Optional normalized version for visual comparison
df1["cci_norm"] = df1["mean_cci"] / (df1["mean_cci"].mean() + EPS)
df2["cci_norm"] = df2["mean_cci"] / (df2["mean_cci"].mean() + EPS)
df3["cci_norm"] = df3["mean_cci"] / (df3["mean_cci"].mean() + EPS)

df1["ratio_norm"] = df1["ratio"] / (df1["ratio"].mean() + EPS)
df2["ratio_norm"] = df2["ratio"] / (df2["ratio"].mean() + EPS)
df3["ratio_norm"] = df3["ratio"] / (df3["ratio"].mean() + EPS)

df1["log_ratio_norm"] = np.log10(df1["ratio_norm"] + EPS)
df2["log_ratio_norm"] = np.log10(df2["ratio_norm"] + EPS)
df3["log_ratio_norm"] = np.log10(df3["ratio_norm"] + EPS)

# ============================================================
# 1. LINEAR COLLAPSE PLOT
# ============================================================

plt.figure(figsize=(7, 5))

plt.scatter(df1["ratio"], df1["mean_cci"], alpha=0.75, label=r"1D ($\alpha=2.5$)")
plt.scatter(df2["ratio"], df2["mean_cci"], alpha=0.75, label=r"2D ($\alpha\gtrsim 3.0$)")
plt.scatter(df3["ratio"], df3["mean_cci"], alpha=0.75, label=r"3D ($\alpha\approx 2.8$ plateau)")

plt.xlabel(r"$R_\alpha = \overline{\dot S}_{cg}^{(+)} / (\overline{I}_{nn}+\varepsilon)^\alpha$")
plt.ylabel("mean CCI")
plt.title("Collapse attempt across dimensions (linear scale)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("collapse_plot_linear.png", dpi=220)
plt.show()

# ============================================================
# 2. LOG-X COLLAPSE PLOT
# ============================================================

plt.figure(figsize=(7, 5))

plt.scatter(df1["log_ratio"], df1["mean_cci"], alpha=0.75, label=r"1D ($\alpha=2.5$)")
plt.scatter(df2["log_ratio"], df2["mean_cci"], alpha=0.75, label=r"2D ($\alpha\gtrsim 3.0$)")
plt.scatter(df3["log_ratio"], df3["mean_cci"], alpha=0.75, label=r"3D ($\alpha\approx 2.8$ plateau)")

plt.xlabel(r"$\log_{10} R_\alpha$")
plt.ylabel("mean CCI")
plt.title("Collapse attempt across dimensions (log-x)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("collapse_plot_logx.png", dpi=220)
plt.show()

# ============================================================
# 3. LOG-LOG COLLAPSE PLOT
# ============================================================

plt.figure(figsize=(7, 5))

plt.scatter(df1["log_ratio"], df1["log_cci"], alpha=0.75, label=r"1D ($\alpha=2.5$)")
plt.scatter(df2["log_ratio"], df2["log_cci"], alpha=0.75, label=r"2D ($\alpha\gtrsim 3.0$)")
plt.scatter(df3["log_ratio"], df3["log_cci"], alpha=0.75, label=r"3D ($\alpha\approx 2.8$ plateau)")

plt.xlabel(r"$\log_{10} R_\alpha$")
plt.ylabel(r"$\log_{10}(\mathrm{CCI})$")
plt.title("Collapse attempt across dimensions (log-log)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("collapse_plot_loglog.png", dpi=220)
plt.show()

# ============================================================
# 4. NORMALIZED COLLAPSE PLOT
# ============================================================

plt.figure(figsize=(7, 5))

plt.scatter(df1["log_ratio_norm"], df1["cci_norm"], alpha=0.75, label=r"1D")
plt.scatter(df2["log_ratio_norm"], df2["cci_norm"], alpha=0.75, label=r"2D")
plt.scatter(df3["log_ratio_norm"], df3["cci_norm"], alpha=0.75, label=r"3D")

plt.xlabel(r"$\log_{10}(R_\alpha / \langle R_\alpha\rangle)$")
plt.ylabel(r"$\mathrm{CCI}/\langle \mathrm{CCI}\rangle$")
plt.title("Normalized collapse across dimensions")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("collapse_plot_normalized.png", dpi=220)
plt.show()

# ============================================================
# 5. OPTIONAL: SIDE-BY-SIDE SUMMARY
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

# left: log-x
ax = axes[0]
ax.scatter(df1["log_ratio"], df1["mean_cci"], alpha=0.75, label="1D")
ax.scatter(df2["log_ratio"], df2["mean_cci"], alpha=0.75, label="2D")
ax.scatter(df3["log_ratio"], df3["mean_cci"], alpha=0.75, label="3D")
ax.set_xlabel(r"$\log_{10} R_\alpha$")
ax.set_ylabel("mean CCI")
ax.set_title("(a) Collapse attempt (log-x)")
ax.grid(alpha=0.3)
ax.legend()

# right: log-log
ax = axes[1]
ax.scatter(df1["log_ratio"], df1["log_cci"], alpha=0.75, label="1D")
ax.scatter(df2["log_ratio"], df2["log_cci"], alpha=0.75, label="2D")
ax.scatter(df3["log_ratio"], df3["log_cci"], alpha=0.75, label="3D")
ax.set_xlabel(r"$\log_{10} R_\alpha$")
ax.set_ylabel(r"$\log_{10}(\mathrm{CCI})$")
ax.set_title("(b) Collapse attempt (log-log)")
ax.grid(alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig("collapse_plot_summary.png", dpi=220)
plt.show()

print("Saved:")
print(" - collapse_plot_linear.png")
print(" - collapse_plot_logx.png")
print(" - collapse_plot_loglog.png")
print(" - collapse_plot_normalized.png")
print(" - collapse_plot_summary.png")
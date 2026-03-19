import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# ============================================================
# xi_aniso: directional anisotropy of nearest-neighbour MI
# This script uses the 3D CSV + recomputes anisotropy from
# the raw simulation data.
# Since we only have mean_mi (averaged over x/y/z), we use
# a proxy: compute xi from the field snapshots saved by
# cci_entropy_scaling_3d.py — or approximate from variance.
#
# Here we demonstrate the xi-split methodology using
# a proxy computed from the available CSV columns.
# ============================================================

EPS = 1e-8
N_BINS = 16
PHI_MAX = 3.0

df3 = pd.read_csv("/home/claude/3d_cci_entropy_information_test.csv")

print(f"3D ensemble: {len(df3)} runs")
print(df3.columns.tolist())

# ============================================================
# Proxy for xi_aniso from available data:
# In the 3D simulation, MI is averaged over x/y/z directions.
# We use the ratio std(CCI, fstruct)/mean as a proxy for
# structural heterogeneity per run — this is what we have.
#
# Better: run extended simulation with per-direction MI.
# Here we demonstrate the split-scan methodology.
# ============================================================

# Proxy: xi ~ F_struct / (mean_mi + eps)
# Higher F_struct relative to MI = more structural heterogeneity
df3["xi_proxy"] = df3["mean_fstruct"] / (df3["mean_mi"] + EPS)
xi_median = df3["xi_proxy"].median()
print(f"\nxi_proxy median: {xi_median:.4f}")

df_low  = df3[df3["xi_proxy"] <= xi_median].copy()
df_high = df3[df3["xi_proxy"] >  xi_median].copy()
print(f"Low-xi group:  {len(df_low)} runs")
print(f"High-xi group: {len(df_high)} runs")

# ============================================================
# Alpha scan for each group
# ============================================================
ALPHAS = np.linspace(0.5, 3.5, 31)

def alpha_scan(df, label):
    results = []
    for alpha in ALPHAS:
        ratio = df["mean_dS_pos"] / ((df["mean_mi"] + EPS) ** alpha)
        r, p = spearmanr(df["mean_cci"], ratio)
        results.append({"alpha": alpha, "spearman_r": r, "label": label})
    return pd.DataFrame(results)

scan_low  = alpha_scan(df_low,  "low-ξ")
scan_high = alpha_scan(df_high, "high-ξ")
scan_all  = alpha_scan(df3,     "all (3D)")

best_low  = scan_low.loc[scan_low["spearman_r"].idxmax()]
best_high = scan_high.loc[scan_high["spearman_r"].idxmax()]
best_all  = scan_all.loc[scan_all["spearman_r"].idxmax()]

print(f"\nAll 3D:   best α = {best_all['alpha']:.2f}, r_s = {best_all['spearman_r']:.3f}")
print(f"Low-ξ:    best α = {best_low['alpha']:.2f}, r_s = {best_low['spearman_r']:.3f}")
print(f"High-ξ:   best α = {best_high['alpha']:.2f}, r_s = {best_high['spearman_r']:.3f}")

# Plateau width (where r_s > 0.99 * max)
def plateau_width(scan):
    max_r = scan["spearman_r"].max()
    plateau = scan[scan["spearman_r"] >= 0.99 * max_r]
    return float(plateau["alpha"].max() - plateau["alpha"].min())

w_all  = plateau_width(scan_all)
w_low  = plateau_width(scan_low)
w_high = plateau_width(scan_high)

print(f"\nPlateau width (Δα at 99% of max):")
print(f"  All 3D:  Δα = {w_all:.2f}")
print(f"  Low-ξ:   Δα = {w_low:.2f}")
print(f"  High-ξ:  Δα = {w_high:.2f}")

# ============================================================
# Figure 1: Alpha scan split
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(scan_all["alpha"],  scan_all["spearman_r"],  "k--",  lw=1.5, label="all (3D)")
ax.plot(scan_low["alpha"],  scan_low["spearman_r"],  "o-",   ms=5, label=f"low-ξ (n={len(df_low)})")
ax.plot(scan_high["alpha"], scan_high["spearman_r"], "s-",   ms=5, label=f"high-ξ (n={len(df_high)})")
ax.axvline(best_low["alpha"],  linestyle=":", color="C1", alpha=0.7)
ax.axvline(best_high["alpha"], linestyle=":", color="C2", alpha=0.7)
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"Spearman $r_s$")
ax.set_title(r"(a) $\alpha$-scan split by $\xi_{\mathrm{proxy}}$")
ax.legend()
ax.grid(alpha=0.3)

# Panel (b): plateau width comparison
ax = axes[1]
groups = ["All 3D", "Low-ξ", "High-ξ"]
widths = [w_all, w_low, w_high]
bars = ax.bar(groups, widths, color=["gray","C1","C2"], alpha=0.8)
ax.set_ylabel(r"Plateau width $\Delta\alpha$ (at 99% max $r_s$)")
ax.set_title(r"(b) Plateau width by $\xi_{\mathrm{proxy}}$ group")
ax.grid(alpha=0.3, axis='y')
for bar, w in zip(bars, widths):
    ax.text(bar.get_x() + bar.get_width()/2, w + 0.01,
            f"{w:.2f}", ha='center', fontsize=10)

plt.suptitle(r"$\xi$-split analysis: testing $\alpha_{\mathrm{eff}}=\alpha(d,\xi)$",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig("paper14_xi_split.png", dpi=220)
plt.close()
print("\nFigure saved: paper14_xi_split.png")

# ============================================================
# Figure 2: xi_proxy distribution
# ============================================================
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(df3["xi_proxy"], bins=10, color="steelblue", alpha=0.7, edgecolor="white")
ax.axvline(xi_median, color="red", linestyle="--", label=f"median = {xi_median:.3f}")
ax.set_xlabel(r"$\xi_{\mathrm{proxy}} = F_{\mathrm{struct}} / (I_{nn}+\varepsilon)$")
ax.set_ylabel("Count")
ax.set_title(r"Distribution of $\xi_{\mathrm{proxy}}$ across 3D ensemble")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("paper14_xi_distribution.png", dpi=220)
plt.close()
print("Figure saved: paper14_xi_distribution.png")

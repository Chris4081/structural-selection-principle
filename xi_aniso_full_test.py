import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression

EPS = 1e-8
ALPHAS = np.linspace(0.5, 3.5, 31)

df = pd.read_csv("/home/claude/3d_cci_entropy_information_v2.csv")
print(f"Loaded v2 CSV: {len(df)} runs")
print(df[["mean_mi_x","mean_mi_y","mean_mi_z","mean_grad_x2","mean_grad_y2","mean_grad_z2"]].describe())

# ============================================================
# Compute xi_aniso (MI-based) and xi_grad (gradient-based)
# ============================================================
mi_xyz = df[["mean_mi_x","mean_mi_y","mean_mi_z"]]
df["xi_aniso"]    = mi_xyz.std(axis=1) / (mi_xyz.mean(axis=1) + EPS)

grad_xyz = df[["mean_grad_x2","mean_grad_y2","mean_grad_z2"]]
df["xi_grad"]     = grad_xyz.std(axis=1) / (grad_xyz.mean(axis=1) + EPS)

# Also old proxy for comparison
df["xi_proxy"]    = df["mean_fstruct"] / (df["mean_mi"] + EPS)

print(f"\nxi_aniso: min={df['xi_aniso'].min():.4f}, max={df['xi_aniso'].max():.4f}, mean={df['xi_aniso'].mean():.4f}")
print(f"xi_grad:  min={df['xi_grad'].min():.4f}, max={df['xi_grad'].max():.4f}, mean={df['xi_grad'].mean():.4f}")

# ============================================================
# Alpha scan per run (for continuous regression)
# ============================================================
def best_alpha_for_run(df_sub):
    """Returns best alpha and plateau width for a subset."""
    results = []
    for alpha in ALPHAS:
        ratio = df_sub["mean_dS_pos"] / ((df_sub["mean_mi"] + EPS) ** alpha)
        if len(df_sub) < 4:
            results.append({"alpha": alpha, "r": 0})
            continue
        r, _ = spearmanr(df_sub["mean_cci"], ratio)
        results.append({"alpha": alpha, "r": r})
    scan = pd.DataFrame(results)
    max_r = scan["r"].max()
    best_a = float(scan.loc[scan["r"].idxmax(), "alpha"])
    plateau = scan[scan["r"] >= 0.99 * max_r]
    width = float(plateau["alpha"].max() - plateau["alpha"].min())
    return best_a, max_r, width

# ============================================================
# Quartil-Split: 4 groups by xi_aniso
# ============================================================
print("\n=== Quartil-Split by xi_aniso ===")
quartiles = pd.qcut(df["xi_aniso"], q=4, labels=["Q1\n(low ξ)", "Q2", "Q3", "Q4\n(high ξ)"])
df["quartile"] = quartiles

quartile_results = []
for q in ["Q1\n(low ξ)", "Q2", "Q3", "Q4\n(high ξ)"]:
    sub = df[df["quartile"] == q]
    alpha, r_s, width = best_alpha_for_run(sub)
    quartile_results.append({"quartile": q, "n": len(sub), "best_alpha": alpha,
                              "r_s": r_s, "plateau_width": width,
                              "mean_xi": sub["xi_aniso"].mean()})
    print(f"  {q.replace(chr(10),' ')}: n={len(sub)}, α={alpha:.2f}, r_s={r_s:.3f}, Δα={width:.2f}, xi_mean={sub['xi_aniso'].mean():.4f}")

qdf = pd.DataFrame(quartile_results)

# ============================================================
# Continuous regression: plateau_width ~ xi_aniso
# ============================================================
print("\n=== Continuous correlation: xi_aniso vs plateau ===")
# For each run, compute a "local" fit quality using leave-one-out alpha
# Simpler: use R² from log-linear fit as proxy for fit quality
results_continuous = []
for _, row in df.iterrows():
    ratio_best = row["mean_dS_pos"] / ((row["mean_mi"] + EPS) ** 2.8)
    results_continuous.append({
        "xi_aniso": row["xi_aniso"],
        "xi_grad": row["xi_grad"],
        "mean_cci": row["mean_cci"],
        "ratio_best": ratio_best,
    })
cont_df = pd.DataFrame(results_continuous)

r_xi_cci, p_xi_cci = spearmanr(cont_df["xi_aniso"], cont_df["mean_cci"])
r_xi_ratio, p_xi_ratio = spearmanr(cont_df["xi_aniso"], cont_df["ratio_best"])
print(f"  Spearman(xi_aniso, CCI): r={r_xi_cci:.3f}, p={p_xi_cci:.4f}")
print(f"  Spearman(xi_aniso, ratio_best): r={r_xi_ratio:.3f}, p={p_xi_ratio:.4f}")

# ============================================================
# Figure 1: Quartil-Split alpha scan
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel (a): alpha scan per quartile
ax = axes[0]
colors_q = ["#2166ac","#4dac26","#d7191c","#a50026"]
labels_q = ["Q1 (low ξ)","Q2","Q3","Q4 (high ξ)"]
for qi, (q_label, color) in enumerate(zip(["Q1\n(low ξ)","Q2","Q3","Q4\n(high ξ)"], colors_q)):
    sub = df[df["quartile"] == q_label]
    scan_q = []
    for alpha in ALPHAS:
        ratio = sub["mean_dS_pos"] / ((sub["mean_mi"] + EPS) ** alpha)
        r, _ = spearmanr(sub["mean_cci"], ratio)
        scan_q.append(r)
    ax.plot(ALPHAS, scan_q, color=color, label=labels_q[qi], linewidth=2)
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"Spearman $r_s$")
ax.set_title(r"(a) $\alpha$-scan by $\xi_{\mathrm{aniso}}$ quartile")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel (b): plateau width by quartile
ax = axes[1]
q_labels_clean = labels_q
widths = [r["plateau_width"] for r in quartile_results]
rs_vals = [r["r_s"] for r in quartile_results]
x = np.arange(4)
bars = ax.bar(x, widths, color=colors_q, alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(q_labels_clean, fontsize=9)
ax.set_ylabel(r"Plateau width $\Delta\alpha$")
ax.set_title(r"(b) Plateau width by $\xi_{\mathrm{aniso}}$ quartile")
ax.grid(alpha=0.3, axis='y')
for bar, w in zip(bars, widths):
    ax.text(bar.get_x()+bar.get_width()/2, w+0.01, f"{w:.2f}", ha='center', fontsize=9)

# Panel (c): continuous scatter xi_aniso vs CCI
ax = axes[2]
ax.scatter(cont_df["xi_aniso"], cont_df["mean_cci"], alpha=0.8, s=60)
z = np.polyfit(cont_df["xi_aniso"], cont_df["mean_cci"], 1)
p = np.poly1d(z)
x_line = np.linspace(cont_df["xi_aniso"].min(), cont_df["xi_aniso"].max(), 50)
ax.plot(x_line, p(x_line), "r--", lw=1.5)
ax.set_xlabel(r"$\xi_{\mathrm{aniso}}$")
ax.set_ylabel("mean CCI")
ax.set_title(fr"(c) CCI vs $\xi_{{\mathrm{{aniso}}}}$" + f"\n(r_s={r_xi_cci:.3f}, p={p_xi_cci:.4f})")
ax.grid(alpha=0.3)

plt.suptitle(r"$\xi_{\mathrm{aniso}}$ analysis: directional MI anisotropy in 3D", fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig("paper14_xi_aniso_quartile.png", dpi=220)
plt.close()
print("Figure saved: paper14_xi_aniso_quartile.png")

# ============================================================
# Figure 2: xi_aniso distribution + comparison
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
ax.hist(df["xi_aniso"], bins=8, color="steelblue", alpha=0.8, edgecolor="white")
for q_val in np.percentile(df["xi_aniso"], [25,50,75]):
    ax.axvline(q_val, color="red", linestyle="--", alpha=0.6)
ax.set_xlabel(r"$\xi_{\mathrm{aniso}} = \mathrm{std}(I_x,I_y,I_z)/\mathrm{mean}(I_x,I_y,I_z)$")
ax.set_ylabel("Count")
ax.set_title(r"Distribution of $\xi_{\mathrm{aniso}}$")
ax.grid(alpha=0.3)

ax = axes[1]
ax.scatter(df["xi_aniso"], df["xi_proxy"]/df["xi_proxy"].max(),
           alpha=0.7, label="xi_proxy (normalized)")
r_comp, _ = spearmanr(df["xi_aniso"], df["xi_proxy"])
ax.set_xlabel(r"$\xi_{\mathrm{aniso}}$ (directional MI)")
ax.set_ylabel(r"$\xi_{\mathrm{proxy}}$ (normalized)")
ax.set_title(fr"Comparison: $r_s={r_comp:.3f}$")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("paper14_xi_comparison.png", dpi=220)
plt.close()
print("Figure saved: paper14_xi_comparison.png")

# Save results
df.to_csv("3d_cci_xi_aniso_results.csv", index=False)
print("Saved: 3d_cci_xi_aniso_results.csv")
print("\nAll done!")

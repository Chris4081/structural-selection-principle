import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib.ticker import ScalarFormatter

# === LOAD DATA ===
df_1d = pd.read_csv("cci_entropy_information_test.csv")
df_2d = pd.read_csv("2d_cci_entropy_information_test.csv")
df_3d = pd.read_csv("3d_cci_entropy_information_test.csv")

datasets = {
    "1D": df_1d,
    "2D": df_2d,
    "3D": df_3d
}

# === ACTUAL COLUMN NAMES FROM YOUR FILES ===
metrics = {
    "CCI": "mean_cci",
    "F_struct": "mean_fstruct",
    "MI": "mean_mi",
    "dS_pos": "mean_dS_pos",
    "ratio": "ratio_best"
}

colors = {
    "1D": "#4C78A8",   # blue
    "2D": "#F58518",   # orange
    "3D": "#54A24B"    # green
}

# === DETECTED FILES ===
print("\n=== DETECTED FILES ===")
for dim, df in datasets.items():
    print(f"{dim}: {df.columns.tolist()}")

# === GLOBAL MEANS ===
print("\n=== GLOBAL MEANS ===")
for metric_name, col in metrics.items():
    for dim, df in datasets.items():
        mean_val = df[col].mean()
        print(f"{dim:>2}  {metric_name:>8}: {mean_val:.6f}")

# === VARIANCE ACROSS DIMENSIONS ===
print("\n=== STABILITY (Variance across dimensions) ===")
for metric_name, col in metrics.items():
    vals = [datasets[dim][col].mean() for dim in ["1D", "2D", "3D"]]
    var = np.var(vals)
    print(f"{metric_name:>8}: variance = {var:.6f}")

# === OPTIONAL INTERNAL CHECK ===
# NOTE: This is not a strong paper metric; keep mainly for internal debugging.
print("\n=== DISTRIBUTION STABILITY (Spearman on sorted values) ===")

def sorted_spearman(a, b):
    a = np.sort(np.asarray(a))
    b = np.sort(np.asarray(b))
    n = min(len(a), len(b))
    return spearmanr(a[:n], b[:n])[0]

for metric_name, col in metrics.items():
    r12 = sorted_spearman(df_1d[col], df_2d[col])
    r23 = sorted_spearman(df_2d[col], df_3d[col])
    r13 = sorted_spearman(df_1d[col], df_3d[col])
    print(f"{metric_name:>8}: 1D-2D = {r12:.3f}, 2D-3D = {r23:.3f}, 1D-3D = {r13:.3f}")

# === PLOT STYLE ===
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
})

# === HISTOGRAM PLOTS ===
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

for i, (metric_name, col) in enumerate(metrics.items()):
    ax = axes[i]

    for dim, df in datasets.items():
        vals = df[col].dropna().values

        ax.hist(
            vals,
            bins=15,
            alpha=0.45,
            color=colors[dim],
            edgecolor="white",
            linewidth=0.8,
            label=dim
        )

        # dashed mean line
        mean_val = np.mean(vals)
        ax.axvline(
            mean_val,
            color=colors[dim],
            linestyle="--",
            linewidth=2,
            alpha=0.95
        )

    ax.set_title(metric_name)
    ax.set_ylabel("Count")
    ax.grid(alpha=0.3, linestyle=":")

    # ratio panel on log scale
    if metric_name == "ratio":
        ax.set_xscale("log")
        ax.set_xlabel("ratio (log scale)")
        ax.xaxis.set_major_formatter(ScalarFormatter())
    else:
        ax.set_xlabel(metric_name)

    ax.legend(frameon=True)

# remove empty subplot if needed
if len(metrics) < len(axes):
    for j in range(len(metrics), len(axes)):
        fig.delaxes(axes[j])

plt.suptitle("Universal Stability Test", fontsize=15, fontweight="bold")
fig.text(
    0.5, 0.01,
    "Dashed lines indicate mean values. The ratio panel is shown on a logarithmic x-scale.",
    ha="center",
    fontsize=10
)

plt.tight_layout(rect=[0, 0.04, 1, 0.95])
plt.savefig("universal_stability_test_final.png", dpi=300, bbox_inches="tight")
plt.savefig("universal_stability_test_final.pdf", bbox_inches="tight")
plt.show()

print("\nSaved plots:")
print(" - universal_stability_test_final.png")
print(" - universal_stability_test_final.pdf")
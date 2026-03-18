import numpy as np
import matplotlib.pyplot as plt

# ============================================
# Current results
# ============================================

dims = np.array([1, 2], dtype=float)

alpha_vals = np.array([2.5, 3.0], dtype=float)
alpha_labels = ["1D: α ≈ 2.5", "2D: α ≳ 3.0"]

spearman_vals = np.array([0.734, 0.870], dtype=float)
spearman_labels = ["1D: r_s = 0.734", "2D: r_s = 0.870"]

# ============================================
# Figure
# ============================================

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# --------------------------------------------
# Panel (a): alpha vs dimension
# --------------------------------------------
ax = axes[0]

# guide to the eye
ax.plot(dims, alpha_vals, linestyle="--", linewidth=1)

# 1D point
ax.scatter([1], [2.5], s=70)
ax.annotate(
    alpha_labels[0],
    xy=(1, 2.5),
    xytext=(1.05, 2.43),
    fontsize=9
)

# 2D lower-bound point
ax.scatter([2], [3.0], s=70, marker="^")
ax.annotate(
    alpha_labels[1],
    xy=(2, 3.0),
    xytext=(1.58, 3.07),
    fontsize=9
)

# upward arrow to indicate lower bound / scan boundary
ax.annotate(
    "",
    xy=(2, 3.18),
    xytext=(2, 3.02),
    arrowprops=dict(arrowstyle="->", lw=1.2)
)

ax.text(1.58, 3.19, "scan boundary", fontsize=8)

ax.set_xlabel("spatial dimension d")
ax.set_ylabel(r"effective exponent $\alpha_{\mathrm{eff}}$")
ax.set_title("(a) Dimensional trend of scaling exponent")
ax.set_xticks([1, 2])
ax.set_ylim(2.3, 3.25)
ax.grid(alpha=0.3)

# --------------------------------------------
# Panel (b): Spearman correlation vs dimension
# --------------------------------------------
ax = axes[1]

# guide to the eye
ax.plot(dims, spearman_vals, linestyle="--", linewidth=1)

ax.scatter(dims, spearman_vals, s=70)

for x, y, lab in zip(dims, spearman_vals, spearman_labels):
    ax.annotate(
        lab,
        xy=(x, y),
        xytext=(x + 0.03, y - 0.03 if x == 2 else y + 0.015),
        fontsize=9
    )

ax.set_xlabel("spatial dimension d")
ax.set_ylabel(r"Spearman correlation $r_s$")
ax.set_title("(b) Robustness of entropy–information scaling")
ax.set_xticks([1, 2])
ax.set_ylim(0.68, 0.90)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("alpha_dimension_figure.png", dpi=200)
plt.show()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA

EPS = 1e-8

# Load all three datasets
df1 = pd.read_csv("/home/claude/cci_entropy_information_test.csv")
df2 = pd.read_csv("/home/claude/2d_cci_entropy_information_test.csv")
df3 = pd.read_csv("/home/claude/3d_cci_entropy_information_test.csv")

# Build log-space coordinates
def logcoords(df):
    return np.column_stack([
        np.log10(df["mean_dS_pos"] + EPS),
        np.log10(df["mean_mi"]     + EPS),
        np.log10(df["mean_cci"]    + EPS),
    ])

X1 = logcoords(df1)
X2 = logcoords(df2)
X3 = logcoords(df3)

colors = {"1D": "#2166ac", "2D": "#d6604d", "3D": "#1a9850"}
markers = {"1D": "o", "2D": "^", "3D": "s"}

# ============================================================
# Figure 1: 3D log-space scatter
# ============================================================
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

for X, label in [(X1,"1D"),(X2,"2D"),(X3,"3D")]:
    ax.scatter(X[:,0], X[:,1], X[:,2],
               c=colors[label], marker=markers[label],
               s=60, alpha=0.8, label=label)

ax.set_xlabel(r"$\log_{10}\,\overline{\dot S}_{cg}^{(+)}$", labelpad=8)
ax.set_ylabel(r"$\log_{10}\,\overline{I}_{nn}$", labelpad=8)
ax.set_zlabel(r"$\log_{10}\,\overline{\mathrm{CCI}}$", labelpad=8)
ax.set_title("Observable space: log-space structure across dimensions",
             pad=14)
ax.legend()
plt.tight_layout()
plt.savefig("paper14_3d_logspace.png", dpi=220)
plt.close()
print("Figure 1 done")

# ============================================================
# Figure 2: PCA — effective dimensionality
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

for idx, (X, label) in enumerate([(X1,"1D"),(X2,"2D"),(X3,"3D")]):
    pca = PCA(n_components=3)
    pca.fit(X)
    var = pca.explained_variance_ratio_
    cum = np.cumsum(var)
    ax = axes[idx]
    bars = ax.bar([1,2,3], var*100, color=colors[label],
                  alpha=0.8, label="individual")
    ax.plot([1,2,3], cum*100, "k--o", ms=5, label="cumulative")
    ax.set_xticks([1,2,3])
    ax.set_xlabel("PCA component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title(f"{label}: PCA decomposition")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # annotate
    for i, v in enumerate(var):
        ax.text(i+1, v*100+2, f"{v*100:.1f}%", ha='center', fontsize=8)

plt.suptitle("PCA of log-space observables across dimensions",
             y=1.02, fontsize=12)
plt.tight_layout()
plt.savefig("paper14_pca.png", dpi=220)
plt.close()
print("Figure 2 done")

# ============================================================
# Figure 3: Combined summary — manifold evidence
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: 2D projection of all dims (PC1 vs PC2 from joint PCA)
ax = axes[0]
X_all = np.vstack([X1, X2, X3])
labels_all = (["1D"]*len(X1) + ["2D"]*len(X2) + ["3D"]*len(X3))
pca_all = PCA(n_components=2)
X_proj = pca_all.fit_transform(X_all)

n1, n2, n3 = len(X1), len(X2), len(X3)
for X_sub, label in [
    (X_proj[:n1], "1D"),
    (X_proj[n1:n1+n2], "2D"),
    (X_proj[n1+n2:], "3D")
]:
    ax.scatter(X_sub[:,0], X_sub[:,1],
               c=colors[label], marker=markers[label],
               s=60, alpha=0.8, label=label)

v = pca_all.explained_variance_ratio_
ax.set_xlabel(f"PC1 ({v[0]*100:.1f}% var)")
ax.set_ylabel(f"PC2 ({v[1]*100:.1f}% var)")
ax.set_title("Joint PCA projection of all dimensions")
ax.legend()
ax.grid(alpha=0.3)

# Right: variance in PC1 per dimension (manifold width)
ax = axes[1]
widths = []
for X in [X1, X2, X3]:
    pca_d = PCA(n_components=3)
    pca_d.fit(X)
    widths.append(pca_d.explained_variance_ratio_)

dims = [1, 2, 3]
w = np.array(widths)
ax.plot(dims, w[:,0]*100, "o-", label="PC1", color="#2166ac", ms=8)
ax.plot(dims, w[:,1]*100, "s--", label="PC2", color="#d6604d", ms=8)
ax.plot(dims, w[:,2]*100, "^:", label="PC3", color="#1a9850", ms=8)
ax.set_xlabel("Spatial dimension d")
ax.set_ylabel("Variance explained (%)")
ax.set_title("PCA spectrum: effective manifold dimensionality")
ax.set_xticks([1,2,3])
ax.legend()
ax.grid(alpha=0.3)

plt.suptitle("Empirical evidence for the scaling manifold",
             y=1.02, fontsize=12)
plt.tight_layout()
plt.savefig("paper14_manifold_evidence.png", dpi=220)
plt.close()
print("Figure 3 done")
print("\nAll figures saved!")

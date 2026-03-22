import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

# =========================================================
# CONFIG
# =========================================================
CSV_PATH = "cci_entropy_information_test.csv"
OUTDIR = "paper17_output"
os.makedirs(OUTDIR, exist_ok=True)

REGIME_ORDER = ["ordered", "critical", "chaotic"]

# ---------------------------------------------------------
# Column aliases adapted to your actual CSVs
# ---------------------------------------------------------
COLUMN_ALIASES = {
    "CCI": [
        "CCI", "CCCI", "cci", "ccci", "mean_cci", "C_CCI"
    ],
    "Fstruct": [
        "Fstruct", "F_struct", "fstruct", "structural_free_energy",
        "mean_fstruct"
    ],
    "Sdot_cg": [
        "Sdot_cg", "Sdot", "Sdotcg", "S_cg_dot", "sdot_cg", "dotS_cg",
        "mean_dS_pos", "dS_pos", "mean_sdot_cg"
    ],
    "Inn": [
        "Inn", "I_nn", "inn", "mutual_information", "Istruct",
        "mean_mi", "mi"
    ],
    "ratio_best": [
        "ratio_best", "best_ratio", "R_alpha", "Ralpha"
    ],
    "regime": [
        "regime", "true_regime", "pred_regime"
    ],
}

# =========================================================
# HELPERS
# =========================================================
def normalize_col(name: str) -> str:
    return name.strip().lower().replace(" ", "").replace("-", "_")

def find_matching_columns(columns):
    norm_map = {normalize_col(c): c for c in columns}
    matches = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        found = None
        for alias in aliases:
            key = normalize_col(alias)
            if key in norm_map:
                found = norm_map[key]
                break
        matches[canonical] = found
    return matches

def classify_threshold(x, t1, t2):
    if x < t1:
        return "ordered"
    elif x < t2:
        return "critical"
    else:
        return "chaotic"

def grid_search_thresholds(values, true_labels, n_grid=200):
    """
    Finds best (t1, t2) for 3-class threshold classifier.
    """
    values = np.asarray(values, dtype=float)
    true_labels = np.asarray(true_labels)

    vmin, vmax = float(np.min(values)), float(np.max(values))
    grid = np.linspace(vmin, vmax, n_grid)

    best = {
        "acc": -1.0,
        "t1": None,
        "t2": None,
        "pred": None,
    }

    for t1 in grid:
        for t2 in grid:
            if t2 <= t1:
                continue
            pred = np.array([classify_threshold(v, t1, t2) for v in values])
            acc = accuracy_score(true_labels, pred)
            if acc > best["acc"]:
                best.update({
                    "acc": acc,
                    "t1": t1,
                    "t2": t2,
                    "pred": pred
                })
    return best

def ordered_cm(y_true, y_pred):
    return confusion_matrix(y_true, y_pred, labels=REGIME_ORDER)

def plot_confusion_matrix(cm, labels, title, savepath):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(cm, interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)

    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center")

    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(savepath, dpi=200, bbox_inches="tight")
    plt.close(fig)

def boxplot_by_regime(df, value_col, regime_col, ylabel, title, savepath):
    data = [df.loc[df[regime_col] == r, value_col].dropna().values for r in REGIME_ORDER]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.boxplot(data, labels=REGIME_ORDER)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(savepath, dpi=200, bbox_inches="tight")
    plt.close(fig)

def scatter_cci_fstruct(df, cci_col, fstruct_col, regime_col, savepath):
    markers = {
        "ordered": "o",
        "critical": "s",
        "chaotic": "^"
    }

    fig, ax = plt.subplots(figsize=(6, 4.8))
    for regime in REGIME_ORDER:
        sub = df[df[regime_col] == regime]
        ax.scatter(
            sub[cci_col],
            sub[fstruct_col],
            label=regime,
            marker=markers[regime]
        )

    ax.set_xlabel("CCI")
    ax.set_ylabel("Structural free energy")
    ax.set_title("CCI vs structural free energy")
    ax.legend()
    fig.tight_layout()
    fig.savefig(savepath, dpi=200, bbox_inches="tight")
    plt.close(fig)

# =========================================================
# LOAD DATA
# =========================================================
df = pd.read_csv(CSV_PATH)
matches = find_matching_columns(df.columns)

cci_col = matches["CCI"]
fstruct_col = matches["Fstruct"]
sdot_col = matches["Sdot_cg"]
inn_col = matches["Inn"]
ratio_col = matches["ratio_best"]
regime_col = matches["regime"]

required = {
    "CCI": cci_col,
    "Fstruct": fstruct_col,
    "Sdot_cg": sdot_col,
    "Inn": inn_col,
    "regime": regime_col,
}

print("=" * 80)
print("Detected columns:")
for k, v in required.items():
    print(f"  {k:10s} -> {v}")

if any(v is None for v in [cci_col, regime_col]):
    raise ValueError("Need at least CCI and regime columns.")

# Keep only relevant rows
df = df.copy()
df = df[df[regime_col].isin(REGIME_ORDER)].copy()
df[cci_col] = pd.to_numeric(df[cci_col], errors="coerce")
df = df.dropna(subset=[cci_col, regime_col])

if fstruct_col is not None:
    df[fstruct_col] = pd.to_numeric(df[fstruct_col], errors="coerce")

# =========================================================
# MAIN ANALYSIS: CCI-only threshold classifier
# =========================================================
best_cci = grid_search_thresholds(df[cci_col].values, df[regime_col].values, n_grid=300)

print("\n" + "=" * 80)
print("BEST CCI-ONLY THRESHOLD CLASSIFIER")
print(f"  tau1 = {best_cci['t1']:.6f}")
print(f"  tau2 = {best_cci['t2']:.6f}")
print(f"  accuracy = {best_cci['acc']:.6f}")

cm_cci = ordered_cm(df[regime_col], best_cci["pred"])
print("\nConfusion matrix (CCI-only):")
print(pd.DataFrame(cm_cci, index=REGIME_ORDER, columns=REGIME_ORDER))

print("\nClassification report (CCI-only):")
print(classification_report(df[regime_col], best_cci["pred"], labels=REGIME_ORDER))

# Save predictions
df["pred_cci"] = best_cci["pred"]

# =========================================================
# OPTIONAL COMPARISON: Fstruct-only threshold classifier
# =========================================================
best_fstruct = None
if fstruct_col is not None and df[fstruct_col].notna().sum() > 0:
    sub = df.dropna(subset=[fstruct_col]).copy()
    best_fstruct = grid_search_thresholds(sub[fstruct_col].values, sub[regime_col].values, n_grid=300)

    print("\n" + "=" * 80)
    print("BEST FSTRUCT-ONLY THRESHOLD CLASSIFIER")
    print(f"  tau1 = {best_fstruct['t1']:.6f}")
    print(f"  tau2 = {best_fstruct['t2']:.6f}")
    print(f"  accuracy = {best_fstruct['acc']:.6f}")

    cm_f = ordered_cm(sub[regime_col], best_fstruct["pred"])
    print("\nConfusion matrix (Fstruct-only):")
    print(pd.DataFrame(cm_f, index=REGIME_ORDER, columns=REGIME_ORDER))

    print("\nClassification report (Fstruct-only):")
    print(classification_report(sub[regime_col], best_fstruct["pred"], labels=REGIME_ORDER))

# =========================================================
# FIGURES
# =========================================================
# Figure 1: CCI by true regime
fig1 = os.path.join(OUTDIR, "paper17_fig1_cci_by_regime.png")
boxplot_by_regime(
    df=df,
    value_col=cci_col,
    regime_col=regime_col,
    ylabel="CCI",
    title="CCI distribution by regime",
    savepath=fig1
)

# Figure 2: Confusion matrix (CCI-only)
fig2 = os.path.join(OUTDIR, "paper17_fig2_confusion_matrix.png")
plot_confusion_matrix(
    cm=cm_cci,
    labels=REGIME_ORDER,
    title="Confusion matrix (CCI-only)",
    savepath=fig2
)

# Figure 3: CCI vs Fstruct
fig3 = None
if fstruct_col is not None:
    fig3 = os.path.join(OUTDIR, "paper17_fig3_cci_fstruct_scatter.png")
    scatter_cci_fstruct(
        df=df.dropna(subset=[fstruct_col]),
        cci_col=cci_col,
        fstruct_col=fstruct_col,
        regime_col=regime_col,
        savepath=fig3
    )

# =========================================================
# SAVE SUMMARY
# =========================================================
summary_path = os.path.join(OUTDIR, "paper17_summary.txt")
with open(summary_path, "w", encoding="utf-8") as f:
    f.write("Paper 17 summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"CSV: {CSV_PATH}\n")
    f.write(f"Rows used: {len(df)}\n\n")

    f.write("Detected columns:\n")
    for k, v in required.items():
        f.write(f"  {k:10s} -> {v}\n")

    f.write("\nCCI-only classifier:\n")
    f.write(f"  tau1 = {best_cci['t1']:.6f}\n")
    f.write(f"  tau2 = {best_cci['t2']:.6f}\n")
    f.write(f"  accuracy = {best_cci['acc']:.6f}\n\n")

    f.write("Confusion matrix (CCI-only):\n")
    f.write(pd.DataFrame(cm_cci, index=REGIME_ORDER, columns=REGIME_ORDER).to_string())
    f.write("\n\n")

    if best_fstruct is not None:
        f.write("Fstruct-only classifier:\n")
        f.write(f"  tau1 = {best_fstruct['t1']:.6f}\n")
        f.write(f"  tau2 = {best_fstruct['t2']:.6f}\n")
        f.write(f"  accuracy = {best_fstruct['acc']:.6f}\n\n")

    f.write("Figures:\n")
    f.write(f"  {fig1}\n")
    f.write(f"  {fig2}\n")
    if fig3 is not None:
        f.write(f"  {fig3}\n")

print("\n" + "=" * 80)
print("Done.")
print(f"Summary written to: {summary_path}")
print("Figures saved in:", OUTDIR)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

EPS = 1e-8

FILES = {
    "1D": "cci_entropy_information_test.csv",
    "2D": "2d_cci_entropy_information_test.csv",
    "3D": "3d_cci_entropy_information_test.csv",
}

results = []

for dim, file in FILES.items():
    df = pd.read_csv(file)

    required_cols = ["mean_cci", "mean_mi", "mean_dS_pos"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"{dim} file missing column: {col}")

    mask = (
        (df["mean_cci"] > 0) &
        (df["mean_mi"] > 0) &
        (df["mean_dS_pos"] > 0)
    )

    X = np.column_stack([
        np.log(df.loc[mask, "mean_dS_pos"].to_numpy()),
        -np.log(df.loc[mask, "mean_mi"].to_numpy())
    ])
    y = np.log(df.loc[mask, "mean_cci"].to_numpy())

    reg = LinearRegression()
    reg.fit(X, y)

    a0 = reg.intercept_
    a = reg.coef_[0]
    b = reg.coef_[1]
    r2 = reg.score(X, y)

    y_pred = reg.predict(X)

    results.append({
        "dimension": dim,
        "intercept": a0,
        "a_dS": a,
        "b_MI": b,
        "R2": r2,
        "n_points": len(y)
    })

    # Save per-dimension fit data
    out_df = pd.DataFrame({
        "log_dS": X[:, 0],
        "minus_log_MI": X[:, 1],
        "log_CCI": y,
        "log_CCI_pred": y_pred
    })
    out_df.to_csv(f"paper13_fit_{dim.lower()}.csv", index=False)

results_df = pd.DataFrame(results)
results_df.to_csv("paper13_multiparameter_results.csv", index=False)

print("\n=== Multi-parameter fit results ===")
print(results_df)

# ------------------------------------------------------------
# Figure 1: coefficients a and b across dimensions
# ------------------------------------------------------------

dims = results_df["dimension"].tolist()
a_vals = results_df["a_dS"].to_numpy()
b_vals = results_df["b_MI"].to_numpy()

x = np.arange(len(dims))

plt.figure(figsize=(7, 5))
plt.plot(x, a_vals, marker="o", label=r"$a$ (entropy production)")
plt.plot(x, b_vals, marker="s", label=r"$b$ (information loss)")
plt.xticks(x, dims)
plt.ylabel("fitted coefficient")
plt.xlabel("dimension")
plt.title("Multi-parameter scaling coefficients across dimensions")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("paper13_ab_coefficients.png", dpi=220)
plt.show()

# ------------------------------------------------------------
# Figure 2: R^2 across dimensions
# ------------------------------------------------------------

r2_vals = results_df["R2"].to_numpy()

plt.figure(figsize=(6, 4.5))
plt.plot(x, r2_vals, marker="o")
plt.xticks(x, dims)
plt.ylabel(r"$R^2$")
plt.xlabel("dimension")
plt.title("Fit quality of multi-parameter scaling model")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("paper13_r2_comparison.png", dpi=220)
plt.show()

# ------------------------------------------------------------
# Figure 3: observed vs predicted log(CCI)
# ------------------------------------------------------------

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=False, sharey=False)

for ax, dim in zip(axes, dims):
    fit_df = pd.read_csv(f"paper13_fit_{dim.lower()}.csv")
    ax.scatter(fit_df["log_CCI_pred"], fit_df["log_CCI"], alpha=0.75)
    lo = min(fit_df["log_CCI_pred"].min(), fit_df["log_CCI"].min())
    hi = max(fit_df["log_CCI_pred"].max(), fit_df["log_CCI"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--")
    ax.set_title(dim)
    ax.set_xlabel(r"predicted $\log(\mathrm{CCI})$")
    ax.set_ylabel(r"observed $\log(\mathrm{CCI})$")
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("paper13_observed_vs_predicted.png", dpi=220)
plt.show()

print("\nSaved:")
print(" - paper13_multiparameter_results.csv")
print(" - paper13_ab_coefficients.png")
print(" - paper13_r2_comparison.png")
print(" - paper13_observed_vs_predicted.png")
print(" - paper13_fit_1d.csv / 2d / 3d")
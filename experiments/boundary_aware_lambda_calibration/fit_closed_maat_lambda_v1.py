# fit_closed_maat_lambda_v1

import json
import numpy as np
from scipy.optimize import minimize

SECTORS = ["H", "B", "S", "V", "R"]

def softmax_neg_energy(defects, lambdas):
    energy = defects @ lambdas
    w = np.exp(-energy - np.max(-energy))
    return w / np.sum(w)

def expected_defects(defects, lambdas):
    p = softmax_neg_energy(defects, lambdas)
    return p @ defects

def loss(log_lambdas, defects, target):
    lambdas = np.exp(log_lambdas)
    exp_d = expected_defects(defects, lambdas)

    fit_loss = np.sum((exp_d - target) ** 2)

    # Regularisierung gegen Kollaps auf nur ein Feld
    lambda_shares = lambdas / np.sum(lambdas)
    uniform = np.ones_like(lambda_shares) / len(lambda_shares)
    balance_loss = 0.02 * np.sum((lambda_shares - uniform) ** 2)

    # hält λ-Werte in realistischer Nähe der bisherigen v12/v13-Kalibrierung
    prior = np.array([1.61087307, 2.42829723, 4.62365496, 4.64491214, 8.83234363])
    prior_loss = 0.001 * np.sum((np.log(lambdas + 1e-12) - np.log(prior)) ** 2)

    return fit_loss + balance_loss + prior_loss

def fit_lambdas(defects, target, init=None):
    defects = np.asarray(defects, dtype=float)
    target = np.asarray(target, dtype=float)

    if init is None:
        init = np.log(np.ones(defects.shape[1]))

    res = minimize(
        loss,
        init,
        args=(defects, target),
        method="L-BFGS-B",
        options={"maxiter": 5000, "ftol": 1e-12}
    )

    lambdas = np.exp(res.x)
    exp_d = expected_defects(defects, lambdas)

    return {
        "success": bool(res.success),
        "loss": float(res.fun),
        "lambdas": dict(zip(SECTORS, lambdas)),
        "shares": dict(zip(SECTORS, lambdas / np.sum(lambdas))),
        "target": dict(zip(SECTORS, target)),
        "expected": dict(zip(SECTORS, exp_d)),
        "lambda_sum": float(np.sum(lambdas)),
    }

# Beispiel: deine v12-Zielwerte aus dem Dokument
target_v12 = np.array([
    0.147084,
    0.140323,
    0.068124,
    0.065946,
    0.030535,
])

# Hier später echte Benchmark-Defekte laden:
# defects = np.loadtxt("maat_defects.csv", delimiter=",", skiprows=1)

if __name__ == "__main__":
    import pandas as pd

    df = pd.read_csv("maat_defects_fused.csv")

    defects = df[["d_H", "d_B", "d_S", "d_V", "d_R"]].to_numpy()

    # Ziel: beste 20% nach Laufzeit = strukturell effizienteste SAT-Instanzen
    if "runtime" in df.columns:
        best = df.nsmallest(int(0.2 * len(df)), "runtime")
        target = best[["d_H", "d_B", "d_S", "d_V", "d_R"]].median().to_numpy()
    else:
        # Fusion-Fallback: beste 20% = niedrigste mittlere Defektlast
        defect_cols = ["d_H", "d_B", "d_S", "d_V", "d_R"]
        df["_defect_mean"] = df[defect_cols].mean(axis=1)
        best = df.nsmallest(int(0.2 * len(df)), "_defect_mean")
        target = best[defect_cols].median().to_numpy()

    result = fit_lambdas(defects, target)

    print("\n=== Closed MAAT Lambda Fit v2 FUSED ===")
    print("success:", result["success"])
    print("loss:", result["loss"])
    print("lambda_sum:", result["lambda_sum"])

    print("\n--- Lambdas ---")
    for k, v in result["lambdas"].items():
        print(f"lambda_{k}: {v:.8f}")

    print("\n--- Shares ---")
    for k, v in result["shares"].items():
        print(f"{k}: {v:.4f}")

    print("\n--- Target vs Expected ---")
    for k in SECTORS:
        print(f"{k}: target={result['target'][k]:.6f}, expected={result['expected'][k]:.6f}")

    outpath = "closed_maat_lambda_fit_results.json"
    with open(outpath, "w") as f:
        json.dump(result, f, indent=2)

    print("\nSaved:", outpath)

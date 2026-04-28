import json
import os
import numpy as np
from scipy.optimize import minimize

# Example sector order
SECTORS = ["H", "B", "S", "V", "R"]
EPS = 1e-12
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def softmax_weights(defects, lambdas):
    """
    defects: shape (N, A), nonnegative defect matrix r_a[X_i]
    lambdas: shape (A,)
    """
    energy = defects @ lambdas
    energy = energy - np.min(energy)  # numerical stability
    w = np.exp(-energy)
    return w / (np.sum(w) + EPS)


def expected_defects(defects, lambdas):
    w = softmax_weights(defects, lambdas)
    return w @ defects


def maxent_loss(log_lambdas, defects, target_means, reg=1e-4):
    """
    Use log_lambdas so lambda_a > 0 automatically.
    """
    lambdas = np.exp(log_lambdas)
    means = expected_defects(defects, lambdas)

    # relative matching loss
    rel = (means - target_means) / (np.abs(target_means) + EPS)
    loss = np.sum(rel ** 2)

    # weak regularisation prevents absurd lambda concentration
    loss += reg * np.sum(log_lambdas ** 2)
    return loss


def calibrate_lambdas(defects, target_means):
    n_sectors = defects.shape[1]

    result = minimize(
        maxent_loss,
        x0=np.zeros(n_sectors),
        args=(defects, target_means),
        method="Nelder-Mead",
        options={
            "maxiter": 20000,
            "xatol": 1e-12,
            "fatol": 1e-12,
        },
    )

    lambdas = np.exp(result.x)
    means = expected_defects(defects, lambdas)
    weights = softmax_weights(defects, lambdas)

    return {
        "success": bool(result.success),
        "loss": float(result.fun),
        "lambdas": lambdas,
        "expected_means": means,
        "weights": weights,
        "optimizer_message": str(result.message),
    }


def make_demo_defects(seed=42, n=500):
    """
    Demo ensemble.
    Replace this later with real defects from your v11 JSON:
    columns = H,B,S,V,R defects.
    """
    rng = np.random.default_rng(seed)

    # Synthetic defect distributions:
    # H/B moderate, S/V often small, R rare large failures
    H = rng.gamma(shape=1.5, scale=0.12, size=n)
    B = rng.gamma(shape=1.3, scale=0.15, size=n)
    S = rng.gamma(shape=1.1, scale=0.08, size=n)
    V = rng.gamma(shape=1.2, scale=0.07, size=n)
    R = rng.gamma(shape=0.9, scale=0.05, size=n)

    # Add a few bad robustness outliers
    outliers = rng.choice(n, size=max(5, n // 20), replace=False)
    R[outliers] += rng.uniform(0.5, 2.0, size=len(outliers))

    return np.vstack([H, B, S, V, R]).T


def main():
    print("\n=== MAAT v12 MaxEnt Lambda Calibration ===\n")

    defects = make_demo_defects()

    # Target means:
    # Option A: use ensemble median defects as neutral MaxEnt targets
    target_means = np.median(defects, axis=0)

    result = calibrate_lambdas(defects, target_means)

    lambdas = result["lambdas"]
    means = result["expected_means"]

    print("Success:", result["success"])
    print("Loss   :", result["loss"])
    print("Message:", result["optimizer_message"])

    print("\n--- Target means vs calibrated expected means ---")
    for name, target, mean in zip(SECTORS, target_means, means):
        print(f"{name}: target={target:.8f}  expected={mean:.8f}  ratio={mean/(target+EPS):.6f}")

    print("\n--- Calibrated lambda_a ---")
    for name, lam in zip(SECTORS, lambdas):
        print(f"lambda_{name}: {lam:.8f}")

    print("\n--- Normalized lambda shares ---")
    shares = lambdas / np.sum(lambdas)
    for name, share in zip(SECTORS, shares):
        print(f"{name}: {share:.4f}")

    output = {
        "sectors": SECTORS,
        "target_means": target_means.tolist(),
        "expected_means": means.tolist(),
        "lambdas": lambdas.tolist(),
        "lambda_shares": shares.tolist(),
        "loss": result["loss"],
        "success": result["success"],
    }

    outpath = os.path.join(BASE_DIR, "naturkonstante_v12_maxent_lambda_results.json")
    with open(outpath, "w") as f:
        json.dump(output, f, indent=2)

    print("\nSaved:", outpath)


if __name__ == "__main__":
    main()

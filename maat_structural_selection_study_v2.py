"""
MAAT Structural Selection Study v2
==================================

New in v2:
- lambda / beta parameter sweep
- preferred activity window instead of monotonic activity reward
- class-level ranking stability analysis
- heatmap-ready outputs
- top-level robustness summary

Requirements:
    pip install numpy pandas matplotlib

Run:
    python maat_structural_selection_study_v2.py
"""

from __future__ import annotations
import json
import os
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Global configuration
# ============================================================

SEED = 42
rng = np.random.default_rng(SEED)

OUTDIR = "maat_results_v2"
os.makedirs(OUTDIR, exist_ok=True)

# Lattice / dynamics
N = 128
DX = 1.0
DT = 0.02
STEPS = 1200
SAVE_EVERY = 10

# Ensemble
RUNS_PER_CLASS = 12
CLASSES = ["vacuum", "kink", "localized", "chaotic"]

# phi^4 potential
def V(phi: np.ndarray) -> np.ndarray:
    return 0.25 * (phi**2 - 1.0) ** 2

def dV(phi: np.ndarray) -> np.ndarray:
    return phi * (phi**2 - 1.0)

# Diagnostic cutoffs
LAMBDA_STAB_CUTOFF = 10.0
LAMBDA_BOUND_CUTOFF = 10.0
PHI_MAX = 3.0

# Activity control
ALPHA_GRAD = 0.5
ACTIVITY_TARGET = 0.42
ACTIVITY_WIDTH = 0.18

# Parameter sweeps
LAMBDA_STAB_LIST = [0.5, 1.0, 2.0]
LAMBDA_CONN_LIST = [0.5, 1.0, 2.0]
LAMBDA_DYN_LIST = [0.5, 1.0, 2.0]
BETA_LIST = [2.0, 4.0, 8.0]

MI_BINS = 16


# ============================================================
# Utilities
# ============================================================

@dataclass
class RunSummary:
    run_id: int
    init_class: str
    ebar: float
    stab_bar: float
    conn_bar: float
    dyn_bar: float
    final_energy_density: float


def laplacian(phi: np.ndarray) -> np.ndarray:
    return (np.roll(phi, -1) + np.roll(phi, 1) - 2.0 * phi) / (DX**2)

def gradient(phi: np.ndarray) -> np.ndarray:
    return (np.roll(phi, -1) - np.roll(phi, 1)) / (2.0 * DX)

def total_energy_density(phi: np.ndarray, pi: np.ndarray) -> float:
    grad = gradient(phi)
    density = 0.5 * pi**2 + 0.5 * grad**2 + V(phi)
    return float(np.mean(density))

def mutual_information_1d(x: np.ndarray, y: np.ndarray, bins: int = 16) -> float:
    hist_xy, _, _ = np.histogram2d(x, y, bins=bins)
    pxy = hist_xy / np.sum(hist_xy)

    px = np.sum(pxy, axis=1, keepdims=True)
    py = np.sum(pxy, axis=0, keepdims=True)

    nz = pxy > 0
    denom = px @ py
    mi = np.sum(pxy[nz] * np.log(pxy[nz] / denom[nz]))
    return float(mi)


# ============================================================
# Initial conditions
# ============================================================

def make_initial_condition(init_class: str) -> Tuple[np.ndarray, np.ndarray]:
    x = np.arange(N) - N / 2

    if init_class == "vacuum":
        sign = rng.choice([-1.0, 1.0])
        phi = sign * np.ones(N) + 0.02 * rng.normal(size=N)
        pi = 0.02 * rng.normal(size=N)

    elif init_class == "kink":
        width = rng.uniform(4.0, 10.0)
        center = rng.uniform(-10, 10)
        phi = np.tanh((x - center) / width)
        if rng.random() < 0.5:
            phi = -phi
        phi += 0.02 * rng.normal(size=N)
        pi = 0.02 * rng.normal(size=N)

    elif init_class == "localized":
        amp = rng.uniform(0.8, 1.8)
        width = rng.uniform(3.0, 8.0)
        center = rng.uniform(-15, 15)
        phi = amp * np.exp(-((x - center) ** 2) / (2.0 * width**2))
        phi += 0.05 * rng.normal(size=N)
        pi = 0.05 * rng.normal(size=N)

    elif init_class == "chaotic":
        phi = rng.normal(0.0, 1.0, size=N)
        pi = rng.normal(0.0, 0.7, size=N)

    else:
        raise ValueError(f"Unknown init_class: {init_class}")

    return phi.astype(float), pi.astype(float)


# ============================================================
# Diagnostics
# ============================================================

def stability_term(phi: np.ndarray, pi: np.ndarray) -> float:
    kin = np.sum(pi**2)
    term_kin = kin / (LAMBDA_STAB_CUTOFF + kin)

    boundary_violation = np.maximum(0.0, np.abs(phi) - PHI_MAX)
    bound = np.sum(boundary_violation**2)
    term_bound = bound / (LAMBDA_BOUND_CUTOFF + bound)

    return float(term_kin + term_bound)

def connectivity_term(phi: np.ndarray, bins: int = 16) -> float:
    x = phi[:-1]
    y = phi[1:]
    mi = mutual_information_1d(x, y, bins=bins)
    return float(-mi)  # high MI lowers E

def activity_score(phi: np.ndarray, pi: np.ndarray) -> float:
    grad = gradient(phi)
    raw = np.sum(pi**2) + ALPHA_GRAD * np.sum(grad**2)
    return float(raw / (10.0 + raw))

def activity_window_penalty(s_dyn: float) -> float:
    return float(((s_dyn - ACTIVITY_TARGET) / ACTIVITY_WIDTH) ** 2)

def structural_energy(
    phi: np.ndarray,
    pi: np.ndarray,
    lambda_stab: float,
    lambda_conn: float,
    lambda_dyn: float,
) -> Tuple[float, float, float, float]:
    e_stab = stability_term(phi, pi)
    e_conn = connectivity_term(phi, bins=MI_BINS)
    s_dyn = activity_score(phi, pi)
    e_dyn = activity_window_penalty(s_dyn)

    e_total = (
        lambda_stab * e_stab
        + lambda_conn * e_conn
        + lambda_dyn * e_dyn
    )

    return float(e_total), float(e_stab), float(e_conn), float(s_dyn)


# ============================================================
# Dynamics
# ============================================================

def acceleration(phi: np.ndarray) -> np.ndarray:
    return laplacian(phi) - dV(phi)

def rk4_step(phi: np.ndarray, pi: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    def f(phi_, pi_):
        return pi_, acceleration(phi_)

    k1_phi, k1_pi = f(phi, pi)
    k2_phi, k2_pi = f(phi + 0.5 * dt * k1_phi, pi + 0.5 * dt * k1_pi)
    k3_phi, k3_pi = f(phi + 0.5 * dt * k2_phi, pi + 0.5 * dt * k2_pi)
    k4_phi, k4_pi = f(phi + dt * k3_phi, pi + dt * k3_pi)

    phi_new = phi + (dt / 6.0) * (k1_phi + 2*k2_phi + 2*k3_phi + k4_phi)
    pi_new = pi + (dt / 6.0) * (k1_pi + 2*k2_pi + 2*k3_pi + k4_pi)
    return phi_new, pi_new


# ============================================================
# Ensemble generation
# ============================================================

def generate_base_ensemble() -> List[Dict]:
    payloads = []
    run_id = 0

    for init_class in CLASSES:
        for _ in range(RUNS_PER_CLASS):
            phi, pi = make_initial_condition(init_class)
            series = []

            for step in range(STEPS):
                phi, pi = rk4_step(phi, pi, DT)
                if step % SAVE_EVERY == 0:
                    series.append({"phi": phi.copy(), "pi": pi.copy()})

            payloads.append({
                "run_id": run_id,
                "init_class": init_class,
                "series": series,
                "phi_final": phi.copy(),
                "pi_final": pi.copy(),
                "final_energy_density": total_energy_density(phi, pi),
            })
            run_id += 1

    return payloads


def score_ensemble(
    payloads: List[Dict],
    lambda_stab: float,
    lambda_conn: float,
    lambda_dyn: float,
    beta: float,
) -> pd.DataFrame:
    rows = []

    for payload in payloads:
        e_list, stab_list, conn_list, dyn_list = [], [], [], []

        for snap in payload["series"]:
            e_total, e_stab, e_conn, s_dyn = structural_energy(
                snap["phi"], snap["pi"],
                lambda_stab=lambda_stab,
                lambda_conn=lambda_conn,
                lambda_dyn=lambda_dyn,
            )
            e_list.append(e_total)
            stab_list.append(e_stab)
            conn_list.append(e_conn)
            dyn_list.append(s_dyn)

        rows.append({
            "run_id": payload["run_id"],
            "init_class": payload["init_class"],
            "ebar": float(np.mean(e_list)),
            "stab_bar": float(np.mean(stab_list)),
            "conn_bar": float(np.mean(conn_list)),
            "dyn_bar": float(np.mean(dyn_list)),
            "final_energy_density": payload["final_energy_density"],
        })

    df = pd.DataFrame(rows)
    raw = np.exp(-beta * df["ebar"].values)
    df["weight"] = raw / np.sum(raw)
    df["lambda_stab"] = lambda_stab
    df["lambda_conn"] = lambda_conn
    df["lambda_dyn"] = lambda_dyn
    df["beta"] = beta
    return df


def class_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("init_class")
        .agg(
            ebar_mean=("ebar", "mean"),
            ebar_std=("ebar", "std"),
            weight_sum=("weight", "sum"),
            stab_mean=("stab_bar", "mean"),
            conn_mean=("conn_bar", "mean"),
            dyn_mean=("dyn_bar", "mean"),
        )
        .reindex(CLASSES)
        .reset_index()
    )


# ============================================================
# Sweep
# ============================================================

def run_sweep(payloads: List[Dict]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    all_runs, all_classes = [], []

    total = len(LAMBDA_STAB_LIST)*len(LAMBDA_CONN_LIST)*len(LAMBDA_DYN_LIST)*len(BETA_LIST)
    done = 0

    for ls in LAMBDA_STAB_LIST:
        for lc in LAMBDA_CONN_LIST:
            for ld in LAMBDA_DYN_LIST:
                for beta in BETA_LIST:
                    df = score_ensemble(payloads, ls, lc, ld, beta)
                    cls = class_summary(df)
                    cls["lambda_stab"] = ls
                    cls["lambda_conn"] = lc
                    cls["lambda_dyn"] = ld
                    cls["beta"] = beta
                    all_runs.append(df)
                    all_classes.append(cls)
                    done += 1
                    if done % 10 == 0:
                        print(f"  Sweep: {done}/{total}")

    return pd.concat(all_runs, ignore_index=True), pd.concat(all_classes, ignore_index=True)


def robustness_table(class_df: pd.DataFrame) -> pd.DataFrame:
    winners = (
        class_df.sort_values(
            ["lambda_stab","lambda_conn","lambda_dyn","beta","weight_sum"],
            ascending=[True,True,True,True,False]
        )
        .groupby(["lambda_stab","lambda_conn","lambda_dyn","beta"])
        .first()
        .reset_index()
    )
    rob = winners["init_class"].value_counts().reindex(CLASSES).fillna(0).reset_index()
    rob.columns = ["init_class", "num_wins"]
    return rob


# ============================================================
# Plots
# ============================================================

def plot_weight_heatmap(class_df, fixed_beta, filename):
    sub = class_df[
        (class_df["init_class"] == "kink") &
        (class_df["lambda_stab"] == 1.0) &
        (class_df["beta"] == fixed_beta)
    ].copy()
    pivot = sub.pivot(index="lambda_dyn", columns="lambda_conn", values="weight_sum")
    plt.figure(figsize=(6,5))
    plt.imshow(pivot.values, origin="lower", aspect="auto")
    plt.xticks(range(len(pivot.columns)), pivot.columns)
    plt.yticks(range(len(pivot.index)), pivot.index)
    plt.xlabel("lambda_conn"); plt.ylabel("lambda_dyn")
    plt.title(f"Kink total weight (beta={fixed_beta}, lambda_stab=1)")
    plt.colorbar(label="weight_sum")
    plt.tight_layout(); plt.savefig(filename, dpi=150); plt.close()

def plot_class_win_bar(robust_df, filename):
    plt.figure(figsize=(7,5))
    plt.bar(robust_df["init_class"], robust_df["num_wins"])
    plt.ylabel("Number of parameter settings won")
    plt.title("Robustness of class dominance across parameter sweep")
    plt.tight_layout(); plt.savefig(filename, dpi=150); plt.close()


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("MAAT Structural Selection Study v2")
    print("=" * 60)
    print(f"Ensemble: {len(CLASSES)} classes × {RUNS_PER_CLASS} runs = {len(CLASSES)*RUNS_PER_CLASS} total")
    print(f"Sweep: {len(LAMBDA_STAB_LIST)}×{len(LAMBDA_CONN_LIST)}×{len(LAMBDA_DYN_LIST)}×{len(BETA_LIST)} = "
          f"{len(LAMBDA_STAB_LIST)*len(LAMBDA_CONN_LIST)*len(LAMBDA_DYN_LIST)*len(BETA_LIST)} settings")
    print()

    print("Step 1: Generating base ensemble...")
    payloads = generate_base_ensemble()
    print(f"  Done: {len(payloads)} trajectories generated.")

    print("Step 2: Parameter sweep...")
    runs_df, class_df = run_sweep(payloads)
    robust_df = robustness_table(class_df)

    runs_df.to_csv(os.path.join(OUTDIR,"runs_sweep.csv"), index=False)
    class_df.to_csv(os.path.join(OUTDIR,"class_sweep.csv"), index=False)
    robust_df.to_csv(os.path.join(OUTDIR,"robustness.csv"), index=False)

    with open(os.path.join(OUTDIR,"ensemble_meta.json"),"w") as f:
        json.dump({
            "seed": SEED, "N": N, "DT": DT, "STEPS": STEPS,
            "runs_per_class": RUNS_PER_CLASS, "classes": CLASSES,
            "lambda_stab_list": LAMBDA_STAB_LIST,
            "lambda_conn_list": LAMBDA_CONN_LIST,
            "lambda_dyn_list": LAMBDA_DYN_LIST,
            "beta_list": BETA_LIST,
            "activity_target": ACTIVITY_TARGET,
            "activity_width": ACTIVITY_WIDTH,
        }, f, indent=2)

    print("Step 3: Generating plots...")
    plot_weight_heatmap(class_df, fixed_beta=4.0,
                        filename=os.path.join(OUTDIR,"heatmap_kink_weight_beta4.png"))
    plot_class_win_bar(robust_df,
                       filename=os.path.join(OUTDIR,"class_dominance.png"))

    print()
    print("=" * 60)
    print("ROBUSTNESS TABLE")
    print("=" * 60)
    print(robust_df.to_string(index=False))

    print()
    print("=" * 60)
    print("CLASS SUMMARY (lambda=1/1/1, beta=4)")
    print("=" * 60)
    ref = class_df[
        (class_df["lambda_stab"]==1.0) &
        (class_df["lambda_conn"]==1.0) &
        (class_df["lambda_dyn"]==1.0) &
        (class_df["beta"]==4.0)
    ]
    print(ref[["init_class","ebar_mean","ebar_std","weight_sum","dyn_mean"]].to_string(index=False))
    print()
    print(f"Results saved to: {OUTDIR}/")

if __name__ == "__main__":
    main()

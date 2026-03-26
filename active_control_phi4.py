"""
MAAT Structural Selection Study v2 - 2D VERSION
s* sweep extension
"""

from __future__ import annotations
import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SEED = 42
rng = np.random.default_rng(SEED)
OUTDIR = "maat_results_2d"
os.makedirs(OUTDIR, exist_ok=True)

N = 48
DX = 1.0
DT = 0.02
STEPS = 800
SAVE_EVERY = 8

RUNS_PER_CLASS = 12
CLASSES = ["vacuum", "domainwall", "localized", "chaotic"]

def V(phi):
    return 0.25 * (phi**2 - 1.0)**2

def dV(phi):
    return phi * (phi**2 - 1.0)

LAMBDA_STAB_CUTOFF = 10.0
LAMBDA_BOUND_CUTOFF = 10.0
PHI_MAX = 3.0
ALPHA_GRAD = 0.5

# Original benchmark
ACTIVITY_TARGET = 0.42
ACTIVITY_WIDTH = 0.18

# New sweep values
ACTIVITY_TARGET_LIST = [0.01, 0.02, 0.05, 0.10, 0.20, 0.42]

LAMBDA_STAB_LIST = [0.5, 1.0, 2.0]
LAMBDA_CONN_LIST = [0.5, 1.0, 2.0]
LAMBDA_DYN_LIST  = [0.5, 1.0, 2.0]
BETA_LIST = [2.0, 4.0, 8.0]
MI_BINS = 16

def laplacian_2d(phi):
    return (
        np.roll(phi, -1, 0)
        + np.roll(phi, 1, 0)
        + np.roll(phi, -1, 1)
        + np.roll(phi, 1, 1)
        - 4 * phi
    ) / (DX**2)

def grad_squared_2d(phi):
    dx = (np.roll(phi, -1, 1) - np.roll(phi, 1, 1)) / (2 * DX)
    dy = (np.roll(phi, -1, 0) - np.roll(phi, 1, 0)) / (2 * DX)
    return dx**2 + dy**2

def total_energy_density(phi, pi):
    return float(np.mean(0.5 * pi**2 + 0.5 * grad_squared_2d(phi) + V(phi)))

def mutual_information_1d(x, y, bins=16):
    hist_xy, _, _ = np.histogram2d(x, y, bins=bins)
    total = np.sum(hist_xy)
    if total <= 0:
        return 0.0

    pxy = hist_xy / total
    px = np.sum(pxy, axis=1, keepdims=True)
    py = np.sum(pxy, axis=0, keepdims=True)
    prod = px @ py
    nz = pxy > 0
    return float(np.sum(pxy[nz] * np.log(pxy[nz] / prod[nz])))

def make_initial_condition(init_class):
    x = np.arange(N) - N / 2
    y = np.arange(N) - N / 2
    X, Y = np.meshgrid(x, y)

    if init_class == "vacuum":
        sign = rng.choice([-1.0, 1.0])
        phi = sign * np.ones((N, N)) + 0.02 * rng.normal(size=(N, N))
        pi = 0.02 * rng.normal(size=(N, N))

    elif init_class == "domainwall":
        w = rng.uniform(4.0, 10.0)
        c = rng.uniform(-10, 10)
        phi = np.tanh((X - c) / w)
        if rng.random() < 0.5:
            phi = -phi
        phi += 0.02 * rng.normal(size=(N, N))
        pi = 0.02 * rng.normal(size=(N, N))

    elif init_class == "localized":
        amp = rng.uniform(0.8, 1.8)
        w = rng.uniform(3.0, 8.0)
        cx, cy = rng.uniform(-15, 15, 2)
        phi = amp * np.exp(-((X - cx)**2 + (Y - cy)**2) / (2 * w**2))
        phi += 0.05 * rng.normal(size=(N, N))
        pi = 0.05 * rng.normal(size=(N, N))

    elif init_class == "chaotic":
        phi = rng.normal(0.0, 1.0, size=(N, N))
        pi = rng.normal(0.0, 0.7, size=(N, N))

    else:
        raise ValueError(f"Unknown class: {init_class}")

    return phi.astype(float), pi.astype(float)

def stability_term(phi, pi):
    kin = np.sum(pi**2)
    t1 = kin / (LAMBDA_STAB_CUTOFF + kin)

    bv = np.maximum(0.0, np.abs(phi) - PHI_MAX)
    b = np.sum(bv**2)
    t2 = b / (LAMBDA_BOUND_CUTOFF + b)

    return float(t1 + t2)

def connectivity_term(phi):
    mi_h = mutual_information_1d(phi[:-1, :].ravel(), phi[1:, :].ravel(), MI_BINS)
    mi_v = mutual_information_1d(phi[:, :-1].ravel(), phi[:, 1:].ravel(), MI_BINS)
    return float(-(mi_h + mi_v) / 2.0)

def activity_score(phi, pi):
    raw = np.mean(pi**2) + ALPHA_GRAD * np.mean(grad_squared_2d(phi))
    return float(raw / (10.0 + raw))

def activity_window_penalty(s, s_target, s_width):
    return float(((s - s_target) / s_width) ** 2)

def structural_energy(phi, pi, ls, lc, ld, s_target, s_width):
    es = stability_term(phi, pi)
    ec = connectivity_term(phi)
    sd = activity_score(phi, pi)
    ed = activity_window_penalty(sd, s_target, s_width)
    return float(ls * es + lc * ec + ld * ed), float(es), float(ec), float(sd)

def rk4_step(phi, pi, dt):
    def accel(p):
        return laplacian_2d(p) - dV(p)

    k1p,  k1pi  = pi, accel(phi)
    k2p,  k2pi  = pi + 0.5 * dt * k1pi, accel(phi + 0.5 * dt * k1p)
    k3p,  k3pi  = pi + 0.5 * dt * k2pi, accel(phi + 0.5 * dt * k2p)
    k4p,  k4pi  = pi + dt * k3pi, accel(phi + dt * k3p)

    phi_new = phi + (dt / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p)
    pi_new  = pi  + (dt / 6.0) * (k1pi + 2 * k2pi + 2 * k3pi + k4pi)
    return phi_new, pi_new

def generate_base_ensemble():
    payloads = []
    run_id = 0

    for cls in CLASSES:
        for _ in range(RUNS_PER_CLASS):
            phi, pi = make_initial_condition(cls)
            series = []

            for step in range(STEPS):
                phi, pi = rk4_step(phi, pi, DT)
                if step % SAVE_EVERY == 0:
                    series.append({"phi": phi.copy(), "pi": pi.copy()})

            payloads.append({
                "run_id": run_id,
                "init_class": cls,
                "series": series,
                "phi_final": phi.copy(),
                "pi_final": pi.copy(),
                "final_energy_density": total_energy_density(phi, pi),
            })
            run_id += 1

    return payloads

def score_ensemble(payloads, ls, lc, ld, beta, s_target, s_width):
    rows = []

    for p in payloads:
        el, sl, cl, dl = [], [], [], []

        for snap in p["series"]:
            e, es, ec, sd = structural_energy(
                snap["phi"], snap["pi"], ls, lc, ld, s_target, s_width
            )
            el.append(e)
            sl.append(es)
            cl.append(ec)
            dl.append(sd)

        rows.append({
            "run_id": p["run_id"],
            "init_class": p["init_class"],
            "ebar": float(np.mean(el)),
            "stab_bar": float(np.mean(sl)),
            "conn_bar": float(np.mean(cl)),
            "dyn_bar": float(np.mean(dl)),
            "final_energy_density": p["final_energy_density"],
        })

    df = pd.DataFrame(rows)
    raw = np.exp(-beta * df["ebar"].values)
    df["weight"] = raw / np.sum(raw)
    df[["lambda_stab", "lambda_conn", "lambda_dyn", "beta", "s_target"]] = (
        ls, lc, ld, beta, s_target
    )
    return df

def class_summary(df):
    return (
        df.groupby("init_class")
        .agg(
            ebar_mean=("ebar", "mean"),
            ebar_std=("ebar", "std"),
            weight_sum=("weight", "sum"),
            dyn_mean=("dyn_bar", "mean"),
        )
        .reindex(CLASSES)
        .reset_index()
    )

def run_s_target_sweep(payloads):
    all_runs = []
    all_classes = []

    total = (
        len(ACTIVITY_TARGET_LIST)
        * len(LAMBDA_STAB_LIST)
        * len(LAMBDA_CONN_LIST)
        * len(LAMBDA_DYN_LIST)
        * len(BETA_LIST)
    )
    done = 0

    for s_target in ACTIVITY_TARGET_LIST:
        for ls in LAMBDA_STAB_LIST:
            for lc in LAMBDA_CONN_LIST:
                for ld in LAMBDA_DYN_LIST:
                    for beta in BETA_LIST:
                        df = score_ensemble(
                            payloads, ls, lc, ld, beta, s_target, ACTIVITY_WIDTH
                        )
                        cls = class_summary(df)
                        cls[["lambda_stab", "lambda_conn", "lambda_dyn", "beta", "s_target"]] = (
                            ls, lc, ld, beta, s_target
                        )

                        all_runs.append(df)
                        all_classes.append(cls)

                        done += 1
                        if done % 25 == 0:
                            print(f"  Sweep progress: {done}/{total}")

    return pd.concat(all_runs, ignore_index=True), pd.concat(all_classes, ignore_index=True)

def robustness_table(class_df):
    winners = (
        class_df.sort_values(
            ["s_target", "lambda_stab", "lambda_conn", "lambda_dyn", "beta", "weight_sum"],
            ascending=[True, True, True, True, True, False],
        )
        .groupby(["s_target", "lambda_stab", "lambda_conn", "lambda_dyn", "beta"])
        .first()
        .reset_index()
    )

    rob = (
        winners.groupby(["s_target", "init_class"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=CLASSES, fill_value=0)
        .reset_index()
    )
    return rob

def plot_weight_vs_s_target(class_df):
    # Focus on the reference lambdas/beta for a clean first plot
    sub = class_df[
        (class_df["lambda_stab"] == 1.0) &
        (class_df["lambda_conn"] == 1.0) &
        (class_df["lambda_dyn"] == 1.0) &
        (class_df["beta"] == 4.0)
    ].copy()

    plt.figure(figsize=(8, 5))

    for cls in CLASSES:
        d = sub[sub["init_class"] == cls].sort_values("s_target")
        plt.plot(d["s_target"], d["weight_sum"], marker="o", label=cls)

    plt.xlabel("Activity target s*")
    plt.ylabel("Class weight sum")
    plt.title("Class weight vs activity target s* (reference setting)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "weight_vs_s_target.png"), dpi=200)
    plt.close()

def plot_winner_phase_diagram(rob_df):
    plt.figure(figsize=(8, 5))

    for cls in CLASSES:
        plt.plot(
            rob_df["s_target"],
            rob_df[cls],
            marker="o",
            label=cls
        )

    plt.xlabel("Activity target s*")
    plt.ylabel("Wins across parameter settings")
    plt.title("Winner counts vs activity target s*")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "winner_counts_vs_s_target.png"), dpi=200)
    plt.close()

def main():
    print("=" * 72)
    print("MAAT STRUCTURAL SELECTION — 2D DOMAIN-WALL STUDY")
    print("with s* sweep")
    print("=" * 72)

    print("Step 1: Generate one fixed ensemble...")
    payloads = generate_base_ensemble()
    print(f"  → {len(payloads)} trajectories generated")

    print("\nStep 2: Sweep over s* and selection parameters...")
    runs_df, class_df = run_s_target_sweep(payloads)

    print("\nStep 3: Build robustness table...")
    rob_df = robustness_table(class_df)

    runs_df.to_csv(os.path.join(OUTDIR, "runs_2d_sstar_sweep.csv"), index=False)
    class_df.to_csv(os.path.join(OUTDIR, "class_2d_sstar_sweep.csv"), index=False)
    rob_df.to_csv(os.path.join(OUTDIR, "robustness_2d_sstar_sweep.csv"), index=False)

    with open(os.path.join(OUTDIR, "meta_2d_sstar_sweep.json"), "w") as f:
        json.dump({
            "N": N,
            "classes": CLASSES,
            "activity_targets": ACTIVITY_TARGET_LIST,
            "activity_width": ACTIVITY_WIDTH,
            "year": 2026,
            "author": "Christof Krieg",
        }, f, indent=2)

    print("\nStep 4: Generate plots...")
    plot_weight_vs_s_target(class_df)
    plot_winner_phase_diagram(rob_df)

    print("\n" + "=" * 72)
    print("REFERENCE SUMMARY (λ=1,1,1 and β=4)")
    print("=" * 72)
    ref = class_df[
        (class_df["lambda_stab"] == 1.0) &
        (class_df["lambda_conn"] == 1.0) &
        (class_df["lambda_dyn"] == 1.0) &
        (class_df["beta"] == 4.0)
    ][["s_target", "init_class", "ebar_mean", "weight_sum", "dyn_mean"]].sort_values(
        ["s_target", "weight_sum"], ascending=[True, False]
    )
    print(ref.round(4).to_string(index=False))

    print(f"\nSaved results to: {OUTDIR}/")

if __name__ == "__main__":
    main()

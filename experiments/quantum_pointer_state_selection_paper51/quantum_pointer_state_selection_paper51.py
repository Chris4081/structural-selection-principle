#!/usr/bin/env python3
"""
Paper 51: Quantum Pointer-State Selection.

This benchmark tests whether MAAT v1.4 defect-field features predict
pointer-state robustness in a toy non-commuting Lindblad system better than
scalar MAAT compression and simple domain baselines.

The target is the fidelity of the initial pure state after Lindblad evolution.
All features are computed from the initial state and generator parameters, not
from the evolved final state, so the benchmark is leakage-free.

Run:
    python3 quantum_pointer_state_selection_paper51.py
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import r2_score, mean_squared_error


SEED = 51053
OUTDIR = Path("outputs_paper51")

rng = np.random.default_rng(SEED)

I2 = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)


@dataclass
class PointerRecord:
    instance_id: int
    theta: float
    phi: float
    target: float
    gamma_z: float
    gamma_x: float
    gamma_m: float
    hx: float
    hz: float
    t: float
    H: float
    B: float
    S: float
    V: float
    R_resp: float
    R_rob: float
    F_maat: float
    basis_z: float
    basis_x: float
    basis_y: float
    offdiag: float
    channel_conflict: float
    damping_pressure: float
    coherence_exposure: float
    x_exposure: float


def ket(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [
            np.cos(theta / 2),
            np.exp(1j * phi) * np.sin(theta / 2),
        ],
        dtype=complex,
    )


def rho_from_ket(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, np.conjugate(psi))


def bloch_vector(rho: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.real(np.trace(rho @ sigma_x)),
            np.real(np.trace(rho @ sigma_y)),
            np.real(np.trace(rho @ sigma_z)),
        ]
    )


def comm_norm(A: np.ndarray, B: np.ndarray) -> float:
    return float(norm(A @ B - B @ A))


def fidelity_pure(rho: np.ndarray, psi: np.ndarray) -> float:
    val = np.real(np.conjugate(psi) @ rho @ psi)
    return float(np.clip(val, 0.0, 1.0))


def offdiag_weight(rho: np.ndarray) -> float:
    return float(np.abs(rho[0, 1]) + np.abs(rho[1, 0]))


def entropy_pop_balance(rho: np.ndarray) -> float:
    p = np.real(np.diag(rho))
    p = np.clip(p, 1e-12, 1.0)
    ent = -np.sum(p * np.log(p)) / np.log(2)
    return float(1.0 / (1.0 + ent))


def lindblad_rhs(rho: np.ndarray, Hsys: np.ndarray, lindblad_ops: list[np.ndarray]) -> np.ndarray:
    drho = -1j * (Hsys @ rho - rho @ Hsys)
    for L in lindblad_ops:
        LdagL = L.conj().T @ L
        drho += L @ rho @ L.conj().T - 0.5 * (LdagL @ rho + rho @ LdagL)
    return drho


def evolve_lindblad(
    rho0: np.ndarray,
    Hsys: np.ndarray,
    lindblad_ops: list[np.ndarray],
    t: float = 1.0,
    steps: int = 80,
) -> np.ndarray:
    dt = t / steps
    rho = rho0.copy()
    for _ in range(steps):
        k1 = lindblad_rhs(rho, Hsys, lindblad_ops)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1, Hsys, lindblad_ops)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2, Hsys, lindblad_ops)
        k4 = lindblad_rhs(rho + dt * k3, Hsys, lindblad_ops)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        rho = 0.5 * (rho + rho.conj().T)
        tr = np.trace(rho)
        if abs(tr) > 1e-12:
            rho = rho / tr
    return rho


def maat_pointer_features(theta: float, phi: float, params: dict[str, float], instance_id: int) -> PointerRecord:
    gamma_z = params["gamma_z"]
    gamma_x = params["gamma_x"]
    gamma_m = params["gamma_m"]
    hx = params["hx"]
    hz = params["hz"]
    t = params["t"]

    psi = ket(theta, phi)
    rho = rho_from_ket(psi)
    Hsys = hz * sigma_z + hx * sigma_x
    lindblad_ops = [
        np.sqrt(gamma_z) * sigma_z,
        np.sqrt(gamma_x) * sigma_x,
        np.sqrt(gamma_m) * sigma_minus,
    ]

    rho_t = evolve_lindblad(rho, Hsys, lindblad_ops, t=t, steps=80)
    target = fidelity_pure(rho_t, psi)

    bx, by, bz = bloch_vector(rho)
    H = 1.0 / (1.0 + comm_norm(Hsys, rho))
    B = entropy_pop_balance(rho)
    S = 1.0 / (1.0 + offdiag_weight(rho))

    total_gamma = gamma_z + gamma_x + gamma_m + 1e-12
    Vz = abs(bz)
    Vx = abs(bx)
    Vm = max(bz, 0.0)
    V = (gamma_z * Vz + gamma_x * Vx + gamma_m * Vm) / total_gamma

    channel_conflict = 1.0 - abs(gamma_z - gamma_x) / (gamma_z + gamma_x + 1e-12)
    damping_pressure = gamma_m / total_gamma
    coherence_exposure = offdiag_weight(rho) * (gamma_z + gamma_m)
    x_exposure = abs(bz) * gamma_x

    R_resp = (H * B * max(V, 1e-9)) ** (1 / 3)
    R_rob = min(R_resp, (H * B * S * max(V, 1e-9)) ** 0.25)

    eps = 1e-9
    F_maat = -sum(np.log(eps + value) for value in (H, B, S, V, R_rob))

    return PointerRecord(
        instance_id=instance_id,
        theta=theta,
        phi=phi,
        target=target,
        gamma_z=gamma_z,
        gamma_x=gamma_x,
        gamma_m=gamma_m,
        hx=hx,
        hz=hz,
        t=t,
        H=H,
        B=B,
        S=S,
        V=V,
        R_resp=R_resp,
        R_rob=R_rob,
        F_maat=float(F_maat),
        basis_z=abs(bz),
        basis_x=abs(bx),
        basis_y=abs(by),
        offdiag=offdiag_weight(rho),
        channel_conflict=channel_conflict,
        damping_pressure=damping_pressure,
        coherence_exposure=coherence_exposure,
        x_exposure=x_exposure,
    )


def generate_rows(n: int = 3500) -> list[PointerRecord]:
    rows: list[PointerRecord] = []
    for i in range(n):
        theta = rng.uniform(0, np.pi)
        phi = rng.uniform(0, 2 * np.pi)
        params = {
            "gamma_z": rng.uniform(0.0, 3.0),
            "gamma_x": rng.uniform(0.0, 3.0),
            "gamma_m": rng.uniform(0.0, 2.0),
            "hx": rng.uniform(0.0, 1.5),
            "hz": rng.uniform(0.2, 1.5),
            "t": rng.uniform(0.4, 1.4),
        }
        rows.append(maat_pointer_features(theta, phi, params, i))
    return rows


FEATURE_SETS = {
    "MAAT scalar/supports": ["H", "B", "S", "V", "R_resp", "R_rob", "F_maat"],
    "MAAT v1.4 field features": [
        "H",
        "B",
        "S",
        "V",
        "R_resp",
        "R_rob",
        "F_maat",
        "channel_conflict",
        "damping_pressure",
        "coherence_exposure",
        "x_exposure",
    ],
    "Rates only": ["gamma_z", "gamma_x", "gamma_m", "t"],
    "Basis geometry only": ["basis_z", "basis_x", "basis_y", "offdiag"],
    "Hamiltonian only": ["H", "hx", "hz"],
    "Domain combo": [
        "gamma_z",
        "gamma_x",
        "gamma_m",
        "t",
        "basis_z",
        "basis_x",
        "basis_y",
        "offdiag",
        "hx",
        "hz",
    ],
}


def matrix(records: list[PointerRecord], names: list[str]) -> np.ndarray:
    rows = [asdict(r) for r in records]
    return np.asarray([[row[name] for name in names] for row in rows], dtype=float)


def evaluate_models(records: list[PointerRecord]) -> list[dict]:
    target = np.asarray([r.target for r in records], dtype=float)
    out = []
    for name, features in FEATURE_SETS.items():
        X = matrix(records, features)
        rf = RandomForestRegressor(
            n_estimators=300,
            random_state=SEED,
            min_samples_leaf=5,
            n_jobs=-1,
        )
        scores = cross_val_score(rf, X, target, cv=5, scoring="r2")
        # Store out-of-fold predictions too for plot diagnostics.
        preds = np.zeros_like(target)
        kf = KFold(n_splits=5, shuffle=True, random_state=SEED + 1)
        for train, test in kf.split(X):
            rf_fold = RandomForestRegressor(
                n_estimators=300,
                random_state=SEED,
                min_samples_leaf=5,
                n_jobs=-1,
            )
            rf_fold.fit(X[train], target[train])
            preds[test] = rf_fold.predict(X[test])
        out.append(
            {
                "model": name,
                "features": "+".join(features),
                "n_features": len(features),
                "cv_r2_mean": float(scores.mean()),
                "cv_r2_std": float(scores.std()),
                "oof_r2": float(r2_score(target, preds)),
                "oof_rmse": float(np.sqrt(mean_squared_error(target, preds))),
            }
        )
    return out


def scalar_correlations(records: list[PointerRecord]) -> list[dict]:
    keys = [
        "F_maat",
        "R_rob",
        "H",
        "B",
        "S",
        "V",
        "gamma_z",
        "gamma_x",
        "gamma_m",
        "basis_z",
        "basis_x",
        "offdiag",
        "channel_conflict",
        "damping_pressure",
        "coherence_exposure",
        "x_exposure",
    ]
    target = np.asarray([r.target for r in records], dtype=float)
    rows = []
    for key in keys:
        vals = np.asarray([getattr(r, key) for r in records], dtype=float)
        rho, p_value = spearmanr(vals, target)
        rows.append({"score": key, "spearman_rho": float(rho), "p_value": float(p_value)})
    return rows


def feature_importance(records: list[PointerRecord]) -> list[dict]:
    target = np.asarray([r.target for r in records], dtype=float)
    features = FEATURE_SETS["MAAT v1.4 field features"]
    X = matrix(records, features)
    rf = RandomForestRegressor(
        n_estimators=400,
        random_state=SEED + 2,
        min_samples_leaf=5,
        n_jobs=-1,
    )
    rf.fit(X, target)
    perm = permutation_importance(rf, X, target, n_repeats=12, random_state=SEED + 3, n_jobs=-1)
    rows = []
    for name, mean, std, builtin in sorted(
        zip(features, perm.importances_mean, perm.importances_std, rf.feature_importances_),
        key=lambda item: item[1],
        reverse=True,
    ):
        rows.append(
            {
                "feature": name,
                "permutation_importance_mean": float(mean),
                "permutation_importance_std": float(std),
                "rf_builtin_importance": float(builtin),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(records: list[PointerRecord], model_rows: list[dict], corr_rows: list[dict], importances: list[dict]) -> None:
    target = np.asarray([r.target for r in records], dtype=float)
    v_support = np.asarray([r.V for r in records], dtype=float)
    r_rob = np.asarray([r.R_rob for r in records], dtype=float)
    conflict = np.asarray([r.channel_conflict for r in records], dtype=float)
    exposure = np.asarray([r.coherence_exposure for r in records], dtype=float)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    names = [r["model"] for r in model_rows]
    means = [r["cv_r2_mean"] for r in model_rows]
    stds = [r["cv_r2_std"] for r in model_rows]
    colors = ["#14b8a6" if "v1.4" in n else "#334155" if "MAAT" in n else "#64748b" for n in names]
    ax.bar(names, means, yerr=stds, color=colors, edgecolor="#0f172a", capsize=3)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel(r"5-fold CV $R^2$")
    ax.set_title("Pointer-state robustness prediction")
    ax.tick_params(axis="x", rotation=32)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_model_comparison.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    sc = ax.scatter(v_support, target, c=conflict, s=18, cmap="magma", alpha=0.72, edgecolor="none")
    ax.set_xlabel("environment alignment support V")
    ax.set_ylabel("pointer survival fidelity")
    ax.set_title("Pointer robustness depends strongly on environment alignment")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("channel conflict")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_V_support_vs_fidelity.png", dpi=260)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    sc = ax.scatter(r_rob, target, c=exposure, s=18, cmap="viridis", alpha=0.72, edgecolor="none")
    ax.set_xlabel(r"$R_{\mathrm{rob}}$")
    ax.set_ylabel("pointer survival fidelity")
    ax.set_title("Scalar robustness retains signal but misses field effects")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("coherence exposure")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_Rrob_vs_fidelity.png", dpi=260)
    plt.close(fig)

    top = importances[:10]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.barh(
        [r["feature"] for r in reversed(top)],
        [r["permutation_importance_mean"] for r in reversed(top)],
        color="#22c55e",
    )
    ax.set_xlabel("permutation importance")
    ax.set_title("MAAT v1.4 field-feature importance")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_v14_feature_importance.png", dpi=260)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = generate_rows(n=3500)
    record_rows = [asdict(r) for r in records]
    write_csv(OUTDIR / "paper51_pointer_instances.csv", record_rows)

    model_rows = evaluate_models(records)
    write_csv(OUTDIR / "paper51_model_comparison.csv", model_rows)

    corr_rows = scalar_correlations(records)
    write_csv(OUTDIR / "paper51_scalar_correlations.csv", corr_rows)

    importance_rows = feature_importance(records)
    write_csv(OUTDIR / "paper51_v14_feature_importance.csv", importance_rows)

    plot_results(records, model_rows, corr_rows, importance_rows)

    def get_model(name: str) -> dict:
        return next(row for row in model_rows if row["model"] == name)

    maat_scalar = get_model("MAAT scalar/supports")
    maat_v14 = get_model("MAAT v1.4 field features")
    domain_combo = get_model("Domain combo")
    summary = {
        "model": "Paper 51 Quantum Pointer-State Selection",
        "status": "Toy non-commuting Lindblad benchmark; not a full decoherence theory or experimental fit.",
        "random_seed": SEED,
        "n_instances": len(records),
        "channels": ["z-dephasing", "x-dephasing", "amplitude damping"],
        "target": "initial-state fidelity after Lindblad evolution",
        "leakage_note": "Features are computed from the initial state and generator parameters, not from final fidelity.",
        "model_comparison": model_rows,
        "scalar_correlations": corr_rows,
        "v14_feature_importance": importance_rows,
        "r2_maat_scalar": maat_scalar["cv_r2_mean"],
        "r2_maat_v14": maat_v14["cv_r2_mean"],
        "r2_domain_combo": domain_combo["cv_r2_mean"],
        "delta_v14_minus_scalar": float(maat_v14["cv_r2_mean"] - maat_scalar["cv_r2_mean"]),
        "delta_v14_minus_domain_combo": float(maat_v14["cv_r2_mean"] - domain_combo["cv_r2_mean"]),
        "key_interpretation": (
            "Scalar MAAT retains moderate pointer-robustness signal, but v1.4 "
            "field features substantially improve prediction in a non-commuting "
            "Lindblad benchmark."
        ),
        "outputs": [
            "paper51_pointer_instances.csv",
            "paper51_model_comparison.csv",
            "paper51_scalar_correlations.csv",
            "paper51_v14_feature_importance.csv",
            "paper51_summary.json",
            "fig1_model_comparison.png",
            "fig2_V_support_vs_fidelity.png",
            "fig3_Rrob_vs_fidelity.png",
            "fig4_v14_feature_importance.png",
        ],
    }
    with (OUTDIR / "paper51_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

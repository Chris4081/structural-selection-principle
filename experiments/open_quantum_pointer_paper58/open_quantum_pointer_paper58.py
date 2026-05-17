#!/usr/bin/env python3
"""
Paper 58: Open quantum pointer stability as a defect-field benchmark.

This is a classical Lindblad simulation. No quantum computer or cloud service
is required. The benchmark asks whether stationarity-sensitive structural
features predict robust pointer-like states better than raw population balance
or simple quantum baselines.

Run:
    python3 open_quantum_pointer_paper58.py
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, LeaveOneGroupOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 58058
N_INSTANCES = 3600
OUTDIR = Path("outputs_paper58")
EPS = 1.0e-10

OUTDIR.mkdir(parents=True, exist_ok=True)
(OUTDIR / ".mplconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUTDIR / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|
SP = SM.conj().T


@dataclass
class QuantumRecord:
    instance_id: int
    family: str
    theta: float
    phi: float
    T: float
    gamma_z: float
    gamma_x: float
    gamma_y: float
    gamma_m: float
    gamma_p: float
    hx: float
    hy: float
    hz: float
    pointer_axis_x: float
    pointer_axis_y: float
    pointer_axis_z: float
    pointer_robustness: float
    fidelity_retention: float
    probability_stability: float
    purity_final: float
    purity_loss: float
    coherence_loss: float
    entropy_final: float
    decoherence_time_frac: float
    H_stationarity: float
    B_pop: float
    B_stat: float
    S_coh: float
    V_align: float
    R_resp_pop: float
    R_rob_pop: float
    F_scalar_pop: float
    R_resp_stat: float
    R_rob_stat: float
    F_scalar_stat: float
    lindblad_norm0: float
    dpdt_norm0: float
    purity_decay_pressure: float
    h_misalignment: float
    channel_conflict: float
    damping_pressure: float
    excitation_pressure: float
    decoherence_pressure: float
    coherence_exposure: float
    transverse_hamiltonian_pressure: float
    pointer_alignment_loss: float
    B_stat_gain: float


def ket(theta: float, phi: float) -> np.ndarray:
    return np.array([np.cos(theta / 2.0), np.exp(1j * phi) * np.sin(theta / 2.0)], dtype=complex)


def rho_from_ket(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def bloch(rho: np.ndarray) -> np.ndarray:
    return np.array(
        [
            float(np.real(np.trace(rho @ SX))),
            float(np.real(np.trace(rho @ SY))),
            float(np.real(np.trace(rho @ SZ))),
        ]
    )


def rho_from_bloch(r: np.ndarray) -> np.ndarray:
    return 0.5 * (I2 + r[0] * SX + r[1] * SY + r[2] * SZ)


def purity(rho: np.ndarray) -> float:
    return float(np.clip(np.real(np.trace(rho @ rho)), 0.0, 1.0))


def entropy_vn(rho: np.ndarray) -> float:
    vals = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
    vals = np.clip(np.real(vals), EPS, 1.0)
    return float(-np.sum(vals * np.log(vals)) / np.log(2.0))


def fidelity_pure(rho: np.ndarray, psi: np.ndarray) -> float:
    return float(np.clip(np.real(psi.conj() @ rho @ psi), 0.0, 1.0))


def lindblad_rhs(rho: np.ndarray, Hsys: np.ndarray, ops: list[np.ndarray]) -> np.ndarray:
    out = -1j * (Hsys @ rho - rho @ Hsys)
    for L in ops:
        LdL = L.conj().T @ L
        out += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return out


def evolve_trajectory(rho0: np.ndarray, Hsys: np.ndarray, ops: list[np.ndarray], T: float, steps: int = 80) -> list[np.ndarray]:
    dt = T / steps
    rho = rho0.copy()
    trajectory = [rho.copy()]
    for _ in range(steps):
        k1 = lindblad_rhs(rho, Hsys, ops)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1, Hsys, ops)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2, Hsys, ops)
        k4 = lindblad_rhs(rho + dt * k3, Hsys, ops)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        rho = 0.5 * (rho + rho.conj().T)
        tr = np.trace(rho)
        if abs(tr) > EPS:
            rho = rho / tr
        # Project tiny numerical eigenvalue violations back to the Bloch ball.
        r = bloch(rho)
        nr = np.linalg.norm(r)
        if nr > 1.0:
            r = r / nr
            rho = rho_from_bloch(r)
        trajectory.append(rho.copy())
    return trajectory


def channel_parameters(family: str, rng: np.random.Generator) -> dict[str, float]:
    # Rates are deliberately overlapping, so the family labels are not trivial.
    if family == "z_dephasing":
        return {
            "gamma_z": rng.uniform(0.35, 2.8),
            "gamma_x": rng.uniform(0.0, 0.25),
            "gamma_y": rng.uniform(0.0, 0.15),
            "gamma_m": rng.uniform(0.0, 0.18),
            "gamma_p": rng.uniform(0.0, 0.08),
        }
    if family == "x_dephasing":
        return {
            "gamma_z": rng.uniform(0.0, 0.22),
            "gamma_x": rng.uniform(0.35, 2.8),
            "gamma_y": rng.uniform(0.0, 0.15),
            "gamma_m": rng.uniform(0.0, 0.14),
            "gamma_p": rng.uniform(0.0, 0.08),
        }
    if family == "amplitude_damping":
        return {
            "gamma_z": rng.uniform(0.0, 0.35),
            "gamma_x": rng.uniform(0.0, 0.20),
            "gamma_y": rng.uniform(0.0, 0.10),
            "gamma_m": rng.uniform(0.45, 2.5),
            "gamma_p": rng.uniform(0.0, 0.08),
        }
    if family == "thermal_relaxation":
        gm = rng.uniform(0.35, 2.2)
        gp = rng.uniform(0.12, 1.2)
        return {
            "gamma_z": rng.uniform(0.0, 0.35),
            "gamma_x": rng.uniform(0.0, 0.22),
            "gamma_y": rng.uniform(0.0, 0.18),
            "gamma_m": gm,
            "gamma_p": gp,
        }
    if family == "mixed_zx":
        base = rng.uniform(0.25, 1.8)
        return {
            "gamma_z": base * rng.uniform(0.65, 1.35),
            "gamma_x": base * rng.uniform(0.65, 1.35),
            "gamma_y": rng.uniform(0.0, 0.25),
            "gamma_m": rng.uniform(0.0, 0.45),
            "gamma_p": rng.uniform(0.0, 0.15),
        }
    if family == "depolarizing":
        base = rng.uniform(0.2, 1.4)
        return {
            "gamma_z": base * rng.uniform(0.75, 1.25),
            "gamma_x": base * rng.uniform(0.75, 1.25),
            "gamma_y": base * rng.uniform(0.75, 1.25),
            "gamma_m": rng.uniform(0.0, 0.18),
            "gamma_p": rng.uniform(0.0, 0.12),
        }
    raise ValueError(f"unknown family {family}")


def pointer_axis(params: dict[str, float]) -> np.ndarray:
    # This is an effective environmental readout axis, not a fundamental claim.
    # Amplitude damping is z-directed, while x/y/z dephasing contributes along
    # the corresponding Pauli axes. Thermal excitation points against damping.
    vec = np.array(
        [
            params["gamma_x"],
            params["gamma_y"],
            params["gamma_z"] + params["gamma_m"] - 0.45 * params["gamma_p"],
        ],
        dtype=float,
    )
    if np.linalg.norm(vec) < EPS:
        return np.array([0.0, 0.0, 1.0])
    return vec / np.linalg.norm(vec)


def make_ops(params: dict[str, float]) -> list[np.ndarray]:
    ops = []
    if params["gamma_z"] > 0:
        ops.append(np.sqrt(params["gamma_z"]) * SZ)
    if params["gamma_x"] > 0:
        ops.append(np.sqrt(params["gamma_x"]) * SX)
    if params["gamma_y"] > 0:
        ops.append(np.sqrt(params["gamma_y"]) * SY)
    if params["gamma_m"] > 0:
        ops.append(np.sqrt(params["gamma_m"]) * SM)
    if params["gamma_p"] > 0:
        ops.append(np.sqrt(params["gamma_p"]) * SP)
    return ops


def compute_record(instance_id: int, family: str, rng: np.random.Generator) -> QuantumRecord:
    theta = float(rng.uniform(0.0, np.pi))
    phi = float(rng.uniform(0.0, 2.0 * np.pi))
    psi = ket(theta, phi)
    rho0 = rho_from_ket(psi)
    params = channel_parameters(family, rng)
    hx = float(rng.uniform(-0.7, 0.7))
    hy = float(rng.uniform(-0.35, 0.35))
    hz = float(rng.uniform(0.15, 1.2))
    T = float(rng.uniform(0.45, 1.6))
    Hsys = hx * SX + hy * SY + hz * SZ
    ops = make_ops(params)
    axis = pointer_axis(params)

    trajectory = evolve_trajectory(rho0, Hsys, ops, T=T, steps=80)
    rhoT = trajectory[-1]
    r0 = bloch(rho0)
    rT = bloch(rhoT)
    rpar0 = float(np.dot(r0, axis))
    rparT = float(np.dot(rT, axis))
    p0 = 0.5 * (1.0 + rpar0)
    pT = 0.5 * (1.0 + rparT)

    rhs0 = lindblad_rhs(rho0, Hsys, ops)
    dr0 = bloch(rhs0)
    dpdt_norm0 = float(abs(0.5 * np.dot(dr0, axis)))
    lindblad_norm0 = float(np.linalg.norm(rhs0, ord="fro"))
    pur0_dot = float(abs(2.0 * np.real(np.trace(rho0 @ rhs0))))

    transverse0 = float(np.sqrt(max(0.0, np.dot(r0, r0) - rpar0**2)))
    denom = math.sqrt(max(EPS, 1.0 - rpar0**2))
    S_coh = float(np.clip(transverse0 / denom, 0.0, 1.0)) if denom > 1e-7 else 0.0
    B_pop = float(1.0 - abs(rpar0))
    B_stat = float(1.0 / (1.0 + 2.0 * dpdt_norm0 + 0.35 * pur0_dot))
    V_align = float(abs(rpar0))
    H_stationarity = float(1.0 / (1.0 + lindblad_norm0))

    total_deph = params["gamma_x"] + params["gamma_y"] + params["gamma_z"]
    total_relax = params["gamma_m"] + params["gamma_p"]
    total_rate = total_deph + total_relax + EPS
    channel_conflict = float(1.0 - max(params["gamma_x"], params["gamma_y"], params["gamma_z"], total_relax) / total_rate)
    damping_pressure = float(params["gamma_m"] / total_rate)
    excitation_pressure = float(params["gamma_p"] / total_rate)
    decoherence_pressure = float(total_rate * T)
    coherence_exposure = float(S_coh * total_rate * T)
    hvec = np.array([hx, hy, hz], dtype=float)
    hnorm = float(np.linalg.norm(hvec))
    transverse_hamiltonian_pressure = float(np.linalg.norm(hvec - np.dot(hvec, axis) * axis) / (hnorm + EPS))
    h_misalignment = transverse_hamiltonian_pressure
    pointer_alignment_loss = float(1.0 - V_align)

    R_resp_pop = float((H_stationarity * max(B_pop, EPS) * max(V_align, EPS)) ** (1.0 / 3.0))
    R_rob_pop = float(min(R_resp_pop, (H_stationarity * max(B_pop, EPS) * max(S_coh, EPS) * max(V_align, EPS)) ** 0.25))
    F_scalar_pop = float(-sum(math.log(EPS + v) for v in [H_stationarity, B_pop, S_coh, V_align, R_rob_pop]))

    R_resp_stat = float((H_stationarity * max(B_stat, EPS) * max(V_align, EPS)) ** (1.0 / 3.0))
    R_rob_stat = float(min(R_resp_stat, (H_stationarity * max(B_stat, EPS) * max(S_coh, EPS) * max(V_align, EPS)) ** 0.25))
    F_scalar_stat = float(-sum(math.log(EPS + v) for v in [H_stationarity, B_stat, S_coh, V_align, R_rob_stat]))
    B_stat_gain = float(B_stat - B_pop)

    fid = fidelity_pure(rhoT, psi)
    prob_stability = float(np.clip(1.0 - abs(pT - p0), 0.0, 1.0))
    purT = purity(rhoT)
    pur_loss = float(1.0 - purT)
    entT = entropy_vn(rhoT)
    transverse_vals = []
    for rho in trajectory:
        r = bloch(rho)
        rp = float(np.dot(r, axis))
        transverse_vals.append(float(np.sqrt(max(0.0, np.dot(r, r) - rp**2))))
    c0 = transverse_vals[0]
    if c0 < 1e-7:
        deco_time = 1.0
    else:
        below = [i for i, c in enumerate(transverse_vals) if c <= 0.5 * c0]
        deco_time = 1.0 if not below else below[0] / max(1, len(transverse_vals) - 1)
    coherence_loss = float(np.clip((c0 - transverse_vals[-1]) / (c0 + EPS), 0.0, 1.0)) if c0 >= 1e-7 else 0.0

    pointer_robustness = float(np.clip(0.50 * fid + 0.35 * prob_stability + 0.15 * purT, 0.0, 1.0))

    return QuantumRecord(
        instance_id=instance_id,
        family=family,
        theta=theta,
        phi=phi,
        T=T,
        gamma_z=params["gamma_z"],
        gamma_x=params["gamma_x"],
        gamma_y=params["gamma_y"],
        gamma_m=params["gamma_m"],
        gamma_p=params["gamma_p"],
        hx=hx,
        hy=hy,
        hz=hz,
        pointer_axis_x=float(axis[0]),
        pointer_axis_y=float(axis[1]),
        pointer_axis_z=float(axis[2]),
        pointer_robustness=pointer_robustness,
        fidelity_retention=fid,
        probability_stability=prob_stability,
        purity_final=purT,
        purity_loss=pur_loss,
        coherence_loss=coherence_loss,
        entropy_final=entT,
        decoherence_time_frac=float(deco_time),
        H_stationarity=H_stationarity,
        B_pop=B_pop,
        B_stat=B_stat,
        S_coh=S_coh,
        V_align=V_align,
        R_resp_pop=R_resp_pop,
        R_rob_pop=R_rob_pop,
        F_scalar_pop=F_scalar_pop,
        R_resp_stat=R_resp_stat,
        R_rob_stat=R_rob_stat,
        F_scalar_stat=F_scalar_stat,
        lindblad_norm0=lindblad_norm0,
        dpdt_norm0=dpdt_norm0,
        purity_decay_pressure=pur0_dot,
        h_misalignment=h_misalignment,
        channel_conflict=channel_conflict,
        damping_pressure=damping_pressure,
        excitation_pressure=excitation_pressure,
        decoherence_pressure=decoherence_pressure,
        coherence_exposure=coherence_exposure,
        transverse_hamiltonian_pressure=transverse_hamiltonian_pressure,
        pointer_alignment_loss=pointer_alignment_loss,
        B_stat_gain=B_stat_gain,
    )


FEATURE_SETS: dict[str, list[str]] = {
    "rates_only": ["gamma_z", "gamma_x", "gamma_y", "gamma_m", "gamma_p", "T"],
    "standard_quantum_baseline": [
        "gamma_z",
        "gamma_x",
        "gamma_y",
        "gamma_m",
        "gamma_p",
        "T",
        "V_align",
        "S_coh",
        "lindblad_norm0",
        "dpdt_norm0",
        "purity_decay_pressure",
        "h_misalignment",
    ],
    "state_geometry_only": ["B_pop", "S_coh", "V_align", "pointer_axis_x", "pointer_axis_y", "pointer_axis_z"],
    "scalar_pop_balance": ["H_stationarity", "B_pop", "S_coh", "V_align", "R_resp_pop", "R_rob_pop", "F_scalar_pop"],
    "scalar_stationarity": ["H_stationarity", "B_stat", "S_coh", "V_align", "R_resp_stat", "R_rob_stat", "F_scalar_stat"],
    "defect_field_stationary": [
        "H_stationarity",
        "B_stat",
        "S_coh",
        "V_align",
        "R_resp_stat",
        "R_rob_stat",
        "F_scalar_stat",
        "lindblad_norm0",
        "dpdt_norm0",
        "purity_decay_pressure",
        "channel_conflict",
        "damping_pressure",
        "excitation_pressure",
        "decoherence_pressure",
        "coherence_exposure",
        "transverse_hamiltonian_pressure",
        "pointer_alignment_loss",
    ],
}


def score(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y, pred).correlation
    return {
        "r2": float(r2_score(y, pred)),
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "spearman": float(0.0 if np.isnan(rho) else rho),
    }


def evaluate_models(df: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = df["pointer_robustness"].to_numpy()
    kfold = KFold(n_splits=5, shuffle=True, random_state=SEED)
    rows = []
    pred_df = df[["instance_id", "family", "pointer_robustness"]].copy()
    for name, cols in FEATURE_SETS.items():
        X = df[cols].to_numpy()
        models = {
            "ridge": make_pipeline(StandardScaler(), Ridge(alpha=2.0)),
            "rf": RandomForestRegressor(n_estimators=140, min_samples_leaf=5, random_state=SEED, n_jobs=1),
        }
        for model_name, model in models.items():
            pred = cross_val_predict(model, X, y, cv=kfold)
            rows.append({"feature_set": name, "model": model_name, "n_features": len(cols), **score(y, pred)})
            if model_name == "rf":
                pred_df[f"pred_{name}"] = pred

    cols = FEATURE_SETS["defect_field_stationary"]
    shuffle_scores = []
    for rep in range(8):
        X = df[cols].to_numpy().copy()
        for j in range(X.shape[1]):
            rng.shuffle(X[:, j])
        model = RandomForestRegressor(n_estimators=100, min_samples_leaf=5, random_state=SEED + rep, n_jobs=1)
        pred = cross_val_predict(model, X, y, cv=kfold)
        shuffle_scores.append(score(y, pred))
    rows.append(
        {
            "feature_set": "shuffled_defect_null",
            "model": "rf",
            "n_features": len(cols),
            "r2": float(np.mean([s["r2"] for s in shuffle_scores])),
            "rmse": float(np.mean([s["rmse"] for s in shuffle_scores])),
            "spearman": float(np.mean([s["spearman"] for s in shuffle_scores])),
        }
    )
    results = pd.DataFrame(rows).sort_values(["model", "r2"], ascending=[True, False])
    return results, pred_df, scalar_correlations(df)


def evaluate_lfo(df: pd.DataFrame) -> pd.DataFrame:
    y = df["pointer_robustness"].to_numpy()
    groups = df["family"].to_numpy()
    logo = LeaveOneGroupOut()
    rows = []
    for name, cols in FEATURE_SETS.items():
        X = df[cols].to_numpy()
        pred = np.full_like(y, fill_value=np.nan, dtype=float)
        for train, test in logo.split(X, y, groups):
            model = RandomForestRegressor(n_estimators=140, min_samples_leaf=5, random_state=SEED, n_jobs=1)
            model.fit(X[train], y[train])
            pred[test] = model.predict(X[test])
            rows.append({"feature_set": name, "held_out_family": groups[test][0], "n_test": int(len(test)), **score(y[test], pred[test])})
        rows.append({"feature_set": name, "held_out_family": "ALL_LFO", "n_test": int(len(y)), **score(y, pred)})
    return pd.DataFrame(rows)


def scalar_correlations(df: pd.DataFrame) -> pd.DataFrame:
    y = df["pointer_robustness"].to_numpy()
    cols = [
        "F_scalar_pop",
        "F_scalar_stat",
        "R_rob_pop",
        "R_rob_stat",
        "B_pop",
        "B_stat",
        "S_coh",
        "V_align",
        "H_stationarity",
        "lindblad_norm0",
        "dpdt_norm0",
        "purity_decay_pressure",
        "coherence_exposure",
        "decoherence_pressure",
        "fidelity_retention",
        "probability_stability",
        "purity_final",
    ]
    rows = []
    for col in cols:
        rho = spearmanr(df[col].to_numpy(), y).correlation
        rows.append({"score": col, "spearman": float(0.0 if np.isnan(rho) else rho), "abs_spearman": float(abs(0.0 if np.isnan(rho) else rho))})
    return pd.DataFrame(rows).sort_values("abs_spearman", ascending=False)


def feature_importance(df: pd.DataFrame) -> pd.DataFrame:
    cols = FEATURE_SETS["defect_field_stationary"]
    X = df[cols].to_numpy()
    y = df["pointer_robustness"].to_numpy()
    model = RandomForestRegressor(n_estimators=180, min_samples_leaf=5, random_state=SEED, n_jobs=1)
    model.fit(X, y)
    perm = permutation_importance(model, X, y, n_repeats=6, random_state=SEED, n_jobs=1)
    return pd.DataFrame(
        {
            "feature": cols,
            "permutation_importance_mean": perm.importances_mean,
            "permutation_importance_std": perm.importances_std,
            "rf_builtin_importance": model.feature_importances_,
        }
    ).sort_values("permutation_importance_mean", ascending=False)


def make_plots(df: pd.DataFrame, results: pd.DataFrame, lfo: pd.DataFrame, corr: pd.DataFrame, importance: pd.DataFrame, pred_df: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    rf = results[results["model"] == "rf"].copy()
    order = rf.sort_values("r2", ascending=False)["feature_set"].tolist()
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.bar(rf.set_index("feature_set").loc[order].index, rf.set_index("feature_set").loc[order, "r2"], color="#2f6f9f")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("5-fold CV R2")
    ax.set_title("Open quantum pointer robustness prediction")
    ax.tick_params(axis="x", rotation=32)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_model_comparison.png", dpi=180)
    plt.close(fig)

    lfo_all = lfo[lfo["held_out_family"] == "ALL_LFO"].copy()
    order = lfo_all.sort_values("spearman", ascending=False)["feature_set"].tolist()
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.bar(lfo_all.set_index("feature_set").loc[order].index, lfo_all.set_index("feature_set").loc[order, "spearman"], color="#ef6f6c")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Leave-channel-family-out Spearman")
    ax.set_title("Channel-family transfer stress test")
    ax.tick_params(axis="x", rotation=32)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_leave_family_out.png", dpi=180)
    plt.close(fig)

    top = importance.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 5.5))
    ax.barh(top["feature"], top["permutation_importance_mean"], xerr=top["permutation_importance_std"], color="#7a5195")
    ax.set_xlabel("Permutation importance")
    ax.set_title("Defect-field stationary feature importance")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_feature_importance.png", dpi=180)
    plt.close(fig)

    topc = corr.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 5.5))
    ax.barh(topc["score"], topc["spearman"], color="#54a24b")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Spearman correlation")
    ax.set_title("Scalar diagnostic correlations")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_scalar_correlations.png", dpi=180)
    plt.close(fig)

    best = rf.sort_values("r2", ascending=False).iloc[0]["feature_set"]
    pred_col = f"pred_{best}"
    fig, ax = plt.subplots(figsize=(6.2, 5.5))
    for fam, sub in pred_df.groupby("family"):
        ax.scatter(sub["pointer_robustness"], sub[pred_col], s=13, alpha=0.65, label=fam)
    lo = min(pred_df["pointer_robustness"].min(), pred_df[pred_col].min())
    hi = max(pred_df["pointer_robustness"].max(), pred_df[pred_col].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1, color="black")
    ax.set_xlabel("Observed pointer robustness")
    ax.set_ylabel(f"Predicted robustness ({best})")
    ax.set_title("Best model prediction scatter")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig5_prediction_scatter.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.scatter(df["B_pop"], df["B_stat"], c=df["pointer_robustness"], cmap="viridis", s=10, alpha=0.7)
    ax.set_xlabel("B_pop = 1 - |p0 - p1|")
    ax.set_ylabel("B_stat = instantaneous stationarity support")
    ax.set_title("Population balance vs stationarity balance")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig6_bpop_vs_bstat.png", dpi=180)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(SEED)
    families = ["z_dephasing", "x_dephasing", "amplitude_damping", "thermal_relaxation", "mixed_zx", "depolarizing"]
    csv_path = OUTDIR / "paper58_pointer_instances.csv"
    if csv_path.exists():
        df_existing = pd.read_csv(csv_path)
        if len(df_existing) >= 1800:
            df = df_existing.sample(n=min(N_INSTANCES, len(df_existing)), random_state=SEED).sort_values("instance_id").reset_index(drop=True)
            print(f"reused {len(df)} existing trajectories from {csv_path}")
        else:
            df = pd.DataFrame()
    else:
        df = pd.DataFrame()
    if df.empty:
        rows = []
        for i in range(N_INSTANCES):
            family = families[i % len(families)]
            rows.append(compute_record(i, family, rng))
            if (i + 1) % 300 == 0:
                print(f"generated {i + 1}/{N_INSTANCES}")
        df = pd.DataFrame([asdict(r) for r in rows])
        df.to_csv(csv_path, index=False)

    results, pred_df, corr = evaluate_models(df, rng)
    lfo = evaluate_lfo(df)
    importance = feature_importance(df)
    results.to_csv(OUTDIR / "paper58_model_comparison.csv", index=False)
    lfo.to_csv(OUTDIR / "paper58_leave_channel_out.csv", index=False)
    corr.to_csv(OUTDIR / "paper58_scalar_correlations.csv", index=False)
    importance.to_csv(OUTDIR / "paper58_feature_importance.csv", index=False)
    pred_df.to_csv(OUTDIR / "paper58_predictions.csv", index=False)
    make_plots(df, results, lfo, corr, importance, pred_df)

    rf = results[results["model"] == "rf"].set_index("feature_set")
    lfo_all = lfo[lfo["held_out_family"] == "ALL_LFO"].set_index("feature_set")
    summary = {
        "model": "Paper 58 Open Quantum Pointer Stability",
        "status": "Classical Lindblad simulation; no quantum computer required; not a collapse theory or hardware result.",
        "seed": SEED,
        "n_instances": int(len(df)),
        "families": families,
        "target": "pointer_robustness = 0.50 fidelity_retention + 0.35 probability_stability + 0.15 purity_final",
        "best_rf_feature_set": str(rf["r2"].idxmax()),
        "rf_results": rf[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "lfo_overall": lfo_all[["r2", "rmse", "spearman"]].to_dict(orient="index"),
        "delta_r2_stationary_vs_population_scalar": float(rf.loc["scalar_stationarity", "r2"] - rf.loc["scalar_pop_balance", "r2"]),
        "delta_r2_defect_vs_standard_quantum": float(rf.loc["defect_field_stationary", "r2"] - rf.loc["standard_quantum_baseline", "r2"]),
        "top_correlations": corr.head(10).to_dict(orient="records"),
        "top_importance": importance.head(10).to_dict(orient="records"),
    }
    with (OUTDIR / "paper58_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

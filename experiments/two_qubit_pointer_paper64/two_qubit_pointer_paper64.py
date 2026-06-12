#!/usr/bin/env python3
"""
Paper 64: Two-qubit pointer selection as an entanglement-robustness benchmark.

This is a fully classical density-matrix Lindblad simulation. No quantum
hardware, cloud backend, or Qiskit runtime is required.

The benchmark asks whether a correlated stationarity balance coordinate
predicts two-qubit entanglement survival better than raw population balance,
and how it compares with simple einselection-inspired baselines.

Run:
    python3 two_qubit_pointer_paper64.py
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

SEED = 64064
N_INSTANCES = 1800
OUTDIR = Path("outputs_paper64")
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

I4 = np.eye(4, dtype=complex)
XX = np.kron(SX, SX)
YY = np.kron(SY, SY)
ZZ = np.kron(SZ, SZ)
ZI = np.kron(SZ, I2)
IZ = np.kron(I2, SZ)
XI = np.kron(SX, I2)
IX = np.kron(I2, SX)
YI = np.kron(SY, I2)
IY = np.kron(I2, SY)
SWAP = np.array(
    [
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
    ],
    dtype=complex,
)
YY2 = np.kron(SY, SY)


@dataclass
class TwoQubitRecord:
    instance_id: int
    state_family: str
    channel_family: str
    T: float
    gamma_lz: float
    gamma_cz: float
    gamma_zz: float
    gamma_amp: float
    gamma_dep: float
    gamma_heat: float
    J: float
    hx1: float
    hx2: float
    hz1: float
    hz2: float
    initial_concurrence: float
    final_concurrence: float
    concurrence_retention: float
    entanglement_robustness: float
    fidelity_retention: float
    purity_final: float
    purity_loss: float
    mutual_info_initial: float
    mutual_info_final: float
    mutual_info_retention: float
    entropy_final: float
    H_stationarity: float
    B_pop_pair: float
    B_corr_stat: float
    S_ent: float
    V_corr: float
    R_resp_pop: float
    R_rob_pop: float
    F_scalar_pop: float
    R_resp_stat: float
    R_rob_stat: float
    F_scalar_stat: float
    lindblad_norm0: float
    local_drift0: float
    corr_drift0: float
    concurrence_slope0: float
    purity_decay0: float
    coherence_norm0: float
    bell_phi_overlap: float
    bell_psi_overlap: float
    dfs_collective_overlap: float
    dfs_parity_overlap: float
    symmetry_expectation: float
    damping_pressure: float
    local_dephasing_pressure: float
    collective_pressure: float
    pair_dephasing_pressure: float
    depolarizing_pressure: float
    hamiltonian_coupling_pressure: float
    stationarity_gain: float


def kron(*ops: np.ndarray) -> np.ndarray:
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def normalize_state(psi: np.ndarray) -> np.ndarray:
    return psi / np.linalg.norm(psi)


def rho_from_ket(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def hermitize_trace_project(rho: np.ndarray) -> np.ndarray:
    rho = 0.5 * (rho + rho.conj().T)
    vals, vecs = np.linalg.eigh(rho)
    vals = np.clip(np.real(vals), 0.0, None)
    if vals.sum() < EPS:
        return I4 / 4.0
    rho = (vecs * vals) @ vecs.conj().T
    return rho / np.trace(rho)


def purity(rho: np.ndarray) -> float:
    return float(np.clip(np.real(np.trace(rho @ rho)), 0.0, 1.0))


def entropy_vn(rho: np.ndarray) -> float:
    vals = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
    vals = np.clip(np.real(vals), EPS, 1.0)
    return float(-np.sum(vals * np.log(vals)) / np.log(2.0))


def partial_trace_2q(rho: np.ndarray, keep: int) -> np.ndarray:
    tensor = rho.reshape(2, 2, 2, 2)
    if keep == 0:
        return np.trace(tensor, axis1=1, axis2=3)
    if keep == 1:
        return np.trace(tensor, axis1=0, axis2=2)
    raise ValueError("keep must be 0 or 1")


def mutual_info(rho: np.ndarray) -> float:
    return max(0.0, entropy_vn(partial_trace_2q(rho, 0)) + entropy_vn(partial_trace_2q(rho, 1)) - entropy_vn(rho))


def concurrence(rho: np.ndarray) -> float:
    rho_h = 0.5 * (rho + rho.conj().T)
    rmat = rho_h @ YY2 @ rho_h.conj() @ YY2
    vals = np.linalg.eigvals(rmat)
    roots = np.sqrt(np.clip(np.real(vals), 0.0, None))
    roots = np.sort(roots)[::-1]
    return float(np.clip(roots[0] - np.sum(roots[1:]), 0.0, 1.0))


def fidelity_pure(rho: np.ndarray, psi: np.ndarray) -> float:
    return float(np.clip(np.real(psi.conj() @ rho @ psi), 0.0, 1.0))


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real(np.trace(rho @ op)))


def offdiag_norm(rho: np.ndarray) -> float:
    off = rho.copy()
    np.fill_diagonal(off, 0.0)
    return float(np.linalg.norm(off, ord="fro"))


def bell_states() -> dict[str, np.ndarray]:
    e00 = np.array([1, 0, 0, 0], dtype=complex)
    e01 = np.array([0, 1, 0, 0], dtype=complex)
    e10 = np.array([0, 0, 1, 0], dtype=complex)
    e11 = np.array([0, 0, 0, 1], dtype=complex)
    return {
        "phi_plus": normalize_state(e00 + e11),
        "phi_minus": normalize_state(e00 - e11),
        "psi_plus": normalize_state(e01 + e10),
        "psi_minus": normalize_state(e01 - e10),
    }


BELL = bell_states()
P_PHI = rho_from_ket(BELL["phi_plus"]) + rho_from_ket(BELL["phi_minus"])
P_PSI = rho_from_ket(BELL["psi_plus"]) + rho_from_ket(BELL["psi_minus"])


def sample_state(family: str, rng: np.random.Generator) -> np.ndarray:
    if family == "product_z":
        idx = int(rng.integers(0, 4))
        psi = np.zeros(4, dtype=complex)
        psi[idx] = 1.0
        # Tiny coherent perturbation prevents degenerate all-zero feature rows.
        psi = normalize_state(0.985 * psi + 0.015 * random_haar_state(rng))
        return psi
    if family == "bell_phi":
        phase = np.exp(1j * rng.uniform(0, 2 * np.pi))
        return normalize_state(np.array([1.0, 0.0, 0.0, phase], dtype=complex))
    if family == "bell_psi":
        phase = np.exp(1j * rng.uniform(0, 2 * np.pi))
        return normalize_state(np.array([0.0, 1.0, phase, 0.0], dtype=complex))
    if family == "partial_phi":
        theta = rng.uniform(0.18 * np.pi, 0.42 * np.pi)
        phase = np.exp(1j * rng.uniform(0, 2 * np.pi))
        return normalize_state(np.array([np.cos(theta), 0.0, 0.0, phase * np.sin(theta)], dtype=complex))
    if family == "partial_psi":
        theta = rng.uniform(0.18 * np.pi, 0.42 * np.pi)
        phase = np.exp(1j * rng.uniform(0, 2 * np.pi))
        return normalize_state(np.array([0.0, np.cos(theta), phase * np.sin(theta), 0.0], dtype=complex))
    if family == "haar_entangled":
        # Reject near-product states so the family tests entanglement survival.
        for _ in range(64):
            psi = random_haar_state(rng)
            if concurrence(rho_from_ket(psi)) > 0.35:
                return psi
        return random_haar_state(rng)
    raise ValueError(f"unknown state family {family}")


def random_haar_state(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=4) + 1j * rng.normal(size=4)
    return normalize_state(z)


def channel_parameters(family: str, rng: np.random.Generator) -> dict[str, float]:
    if family == "local_z_dephasing":
        return {
            "gamma_lz": rng.uniform(0.45, 2.2),
            "gamma_cz": rng.uniform(0.0, 0.25),
            "gamma_zz": rng.uniform(0.0, 0.12),
            "gamma_amp": rng.uniform(0.0, 0.10),
            "gamma_dep": rng.uniform(0.0, 0.08),
            "gamma_heat": rng.uniform(0.0, 0.04),
            "J": rng.uniform(0.0, 0.20),
        }
    if family == "collective_z_dephasing":
        return {
            "gamma_lz": rng.uniform(0.0, 0.18),
            "gamma_cz": rng.uniform(0.45, 2.2),
            "gamma_zz": rng.uniform(0.0, 0.15),
            "gamma_amp": rng.uniform(0.0, 0.08),
            "gamma_dep": rng.uniform(0.0, 0.08),
            "gamma_heat": rng.uniform(0.0, 0.04),
            "J": rng.uniform(0.0, 0.25),
        }
    if family == "pair_zz_dephasing":
        return {
            "gamma_lz": rng.uniform(0.0, 0.16),
            "gamma_cz": rng.uniform(0.0, 0.24),
            "gamma_zz": rng.uniform(0.50, 2.5),
            "gamma_amp": rng.uniform(0.0, 0.08),
            "gamma_dep": rng.uniform(0.0, 0.08),
            "gamma_heat": rng.uniform(0.0, 0.04),
            "J": rng.uniform(0.0, 0.22),
        }
    if family == "local_amplitude_damping":
        return {
            "gamma_lz": rng.uniform(0.0, 0.25),
            "gamma_cz": rng.uniform(0.0, 0.15),
            "gamma_zz": rng.uniform(0.0, 0.10),
            "gamma_amp": rng.uniform(0.45, 2.4),
            "gamma_dep": rng.uniform(0.0, 0.10),
            "gamma_heat": rng.uniform(0.0, 0.08),
            "J": rng.uniform(0.0, 0.18),
        }
    if family == "thermal_mixed":
        return {
            "gamma_lz": rng.uniform(0.0, 0.35),
            "gamma_cz": rng.uniform(0.0, 0.25),
            "gamma_zz": rng.uniform(0.0, 0.16),
            "gamma_amp": rng.uniform(0.20, 1.4),
            "gamma_dep": rng.uniform(0.0, 0.16),
            "gamma_heat": rng.uniform(0.12, 1.0),
            "J": rng.uniform(0.0, 0.25),
        }
    if family == "depolarizing":
        base = rng.uniform(0.25, 1.4)
        return {
            "gamma_lz": rng.uniform(0.0, 0.18),
            "gamma_cz": rng.uniform(0.0, 0.18),
            "gamma_zz": rng.uniform(0.0, 0.12),
            "gamma_amp": rng.uniform(0.0, 0.12),
            "gamma_dep": base,
            "gamma_heat": rng.uniform(0.0, 0.08),
            "J": rng.uniform(0.0, 0.20),
        }
    raise ValueError(f"unknown channel family {family}")


def make_hamiltonian(params: dict[str, float], rng: np.random.Generator) -> tuple[np.ndarray, dict[str, float]]:
    hx1 = rng.uniform(-0.22, 0.22)
    hx2 = rng.uniform(-0.22, 0.22)
    hz1 = rng.uniform(-0.18, 0.18)
    hz2 = rng.uniform(-0.18, 0.18)
    j = params["J"]
    H = hx1 * XI + hx2 * IX + hz1 * ZI + hz2 * IZ + j * (XX + YY)
    return H, {"hx1": hx1, "hx2": hx2, "hz1": hz1, "hz2": hz2}


def make_ops(params: dict[str, float]) -> list[np.ndarray]:
    ops: list[np.ndarray] = []
    if params["gamma_lz"] > 0:
        g = math.sqrt(params["gamma_lz"])
        ops.append(g * ZI)
        ops.append(g * IZ)
    if params["gamma_cz"] > 0:
        ops.append(math.sqrt(params["gamma_cz"]) * (ZI + IZ) / math.sqrt(2.0))
    if params["gamma_zz"] > 0:
        ops.append(math.sqrt(params["gamma_zz"]) * ZZ)
    if params["gamma_amp"] > 0:
        g = math.sqrt(params["gamma_amp"])
        ops.append(g * kron(SM, I2))
        ops.append(g * kron(I2, SM))
    if params["gamma_heat"] > 0:
        g = math.sqrt(params["gamma_heat"])
        ops.append(g * kron(SP, I2))
        ops.append(g * kron(I2, SP))
    if params["gamma_dep"] > 0:
        g = math.sqrt(params["gamma_dep"] / 3.0)
        for op in [XI, IX, YI, IY, ZI, IZ]:
            ops.append(g * op)
    return ops


def lindblad_rhs(rho: np.ndarray, Hsys: np.ndarray, ops: list[np.ndarray]) -> np.ndarray:
    out = -1j * (Hsys @ rho - rho @ Hsys)
    for L in ops:
        LdL = L.conj().T @ L
        out += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return out


def evolve_trajectory(rho0: np.ndarray, Hsys: np.ndarray, ops: list[np.ndarray], T: float, steps: int = 70) -> list[np.ndarray]:
    dt = T / steps
    rho = rho0.copy()
    trajectory = [rho.copy()]
    for _ in range(steps):
        k1 = lindblad_rhs(rho, Hsys, ops)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1, Hsys, ops)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2, Hsys, ops)
        k4 = lindblad_rhs(rho + dt * k3, Hsys, ops)
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        rho = hermitize_trace_project(rho)
        trajectory.append(rho.copy())
    return trajectory


def support(x: float, scale: float = 1.0) -> float:
    return float(1.0 / (1.0 + max(0.0, x) / max(scale, EPS)))


def overlap(rho: np.ndarray, projector: np.ndarray) -> float:
    return float(np.clip(np.real(np.trace(rho @ projector)), 0.0, 1.0))


def generator_features(rho0: np.ndarray, Hsys: np.ndarray, ops: list[np.ndarray], params: dict[str, float]) -> dict[str, float]:
    Lrho = lindblad_rhs(rho0, Hsys, ops)
    rho_eps = hermitize_trace_project(rho0 + 1.0e-3 * Lrho)
    local_drift = abs(np.real(np.trace(Lrho @ ZI))) + abs(np.real(np.trace(Lrho @ IZ)))
    corr_ops = [ZZ, XX, YY]
    corr_drift = sum(abs(np.real(np.trace(Lrho @ op))) for op in corr_ops) / len(corr_ops)
    c0 = concurrence(rho0)
    c_eps = concurrence(rho_eps)
    c_slope = max(0.0, c0 - c_eps) / 1.0e-3
    purity_decay = max(0.0, purity(rho0) - purity(rho_eps)) / 1.0e-3
    lindblad_norm = float(np.linalg.norm(Lrho, ord="fro"))
    coh = offdiag_norm(rho0)
    phi_ov = overlap(rho0, P_PHI)
    psi_ov = overlap(rho0, P_PSI)
    dfs_collective = psi_ov
    dfs_parity = max(phi_ov, psi_ov)
    symmetry = float(np.clip(np.real(np.trace(rho0 @ SWAP)), -1.0, 1.0))
    local_deph = params["gamma_lz"] * (phi_ov + psi_ov)
    collective_protection = params["gamma_cz"] * dfs_collective
    pair_protection = params["gamma_zz"] * dfs_parity
    damping_pressure = params["gamma_amp"] * (1.0 + abs(expect(rho0, ZI) + expect(rho0, IZ)) / 2.0)
    depol_pressure = params["gamma_dep"] * max(0.05, c0)
    return {
        "lindblad_norm0": lindblad_norm,
        "local_drift0": float(local_drift),
        "corr_drift0": float(corr_drift),
        "concurrence_slope0": float(c_slope),
        "purity_decay0": float(purity_decay),
        "coherence_norm0": float(coh),
        "bell_phi_overlap": float(phi_ov),
        "bell_psi_overlap": float(psi_ov),
        "dfs_collective_overlap": float(dfs_collective),
        "dfs_parity_overlap": float(dfs_parity),
        "symmetry_expectation": float(symmetry),
        "damping_pressure": float(damping_pressure),
        "local_dephasing_pressure": float(local_deph),
        "collective_pressure": float(collective_protection),
        "pair_dephasing_pressure": float(pair_protection),
        "depolarizing_pressure": float(depol_pressure),
        "hamiltonian_coupling_pressure": float(abs(params["J"]) * (1.0 + c0)),
    }


def make_record(instance_id: int, rng: np.random.Generator) -> TwoQubitRecord:
    state_families = ["product_z", "bell_phi", "bell_psi", "partial_phi", "partial_psi", "haar_entangled"]
    channel_families = [
        "local_z_dephasing",
        "collective_z_dephasing",
        "pair_zz_dephasing",
        "local_amplitude_damping",
        "thermal_mixed",
        "depolarizing",
    ]
    state_family = state_families[instance_id % len(state_families)]
    channel_family = channel_families[(instance_id // len(state_families)) % len(channel_families)]
    # Add mild randomization after stratified cycling.
    if rng.random() < 0.35:
        state_family = str(rng.choice(state_families))
    if rng.random() < 0.35:
        channel_family = str(rng.choice(channel_families))

    psi = sample_state(state_family, rng)
    rho0 = rho_from_ket(psi)
    params = channel_parameters(channel_family, rng)
    Hsys, hpars = make_hamiltonian(params, rng)
    ops = make_ops(params)
    T = float(rng.uniform(0.55, 1.75))
    traj = evolve_trajectory(rho0, Hsys, ops, T=T)
    rhoT = traj[-1]

    c0 = concurrence(rho0)
    cT = concurrence(rhoT)
    retention = float(np.clip(cT / max(c0, 0.08), 0.0, 1.0))
    fid = fidelity_pure(rhoT, psi)
    pf = purity(rhoT)
    mi0 = mutual_info(rho0)
    miT = mutual_info(rhoT)
    mi_ret = float(np.clip(miT / max(mi0, 0.08), 0.0, 1.0))
    y = float(np.clip(0.58 * cT + 0.22 * retention + 0.10 * mi_ret + 0.10 * pf, 0.0, 1.0))

    feats = generator_features(rho0, Hsys, ops, params)
    H_stat = support(feats["lindblad_norm0"], scale=1.8)
    p1 = (1.0 + expect(rho0, ZI)) / 2.0
    p2 = (1.0 + expect(rho0, IZ)) / 2.0
    B_pop = float(np.clip((1.0 - abs(2.0 * p1 - 1.0)) * (1.0 - abs(2.0 * p2 - 1.0)), 0.0, 1.0))
    B_corr_stat = support(
        0.65 * feats["local_drift0"]
        + 1.15 * feats["corr_drift0"]
        + 0.90 * feats["concurrence_slope0"]
        + 0.25 * feats["purity_decay0"],
        scale=2.2,
    )
    S_ent = float(np.clip(c0, 0.0, 1.0))
    protection = 0.25 * feats["bell_phi_overlap"] + 0.25 * feats["bell_psi_overlap"]
    protection += 0.30 * min(1.0, feats["collective_pressure"]) * feats["dfs_collective_overlap"]
    protection += 0.30 * min(1.0, feats["pair_dephasing_pressure"]) * feats["dfs_parity_overlap"]
    fragility = feats["local_dephasing_pressure"] + feats["damping_pressure"] + feats["depolarizing_pressure"]
    V_corr = float(np.clip((protection + 0.15 * mi0) / (1.0 + 0.45 * fragility), 0.0, 1.0))

    R_resp_pop = float((max(H_stat, EPS) * max(B_pop, EPS) * max(V_corr, EPS)) ** (1.0 / 3.0))
    R_rob_pop = float(min(R_resp_pop, (max(H_stat, EPS) * max(B_pop, EPS) * max(S_ent, EPS) * max(V_corr, EPS)) ** 0.25))
    F_pop = float(-np.log(EPS + H_stat) - np.log(EPS + B_pop) - np.log(EPS + S_ent) - np.log(EPS + V_corr) - np.log(EPS + R_rob_pop))

    R_resp_stat = float((max(H_stat, EPS) * max(B_corr_stat, EPS) * max(V_corr, EPS)) ** (1.0 / 3.0))
    R_rob_stat = float(min(R_resp_stat, (max(H_stat, EPS) * max(B_corr_stat, EPS) * max(S_ent, EPS) * max(V_corr, EPS)) ** 0.25))
    F_stat = float(-np.log(EPS + H_stat) - np.log(EPS + B_corr_stat) - np.log(EPS + S_ent) - np.log(EPS + V_corr) - np.log(EPS + R_rob_stat))

    return TwoQubitRecord(
        instance_id=instance_id,
        state_family=state_family,
        channel_family=channel_family,
        T=T,
        gamma_lz=params["gamma_lz"],
        gamma_cz=params["gamma_cz"],
        gamma_zz=params["gamma_zz"],
        gamma_amp=params["gamma_amp"],
        gamma_dep=params["gamma_dep"],
        gamma_heat=params["gamma_heat"],
        J=params["J"],
        hx1=hpars["hx1"],
        hx2=hpars["hx2"],
        hz1=hpars["hz1"],
        hz2=hpars["hz2"],
        initial_concurrence=c0,
        final_concurrence=cT,
        concurrence_retention=retention,
        entanglement_robustness=y,
        fidelity_retention=fid,
        purity_final=pf,
        purity_loss=max(0.0, 1.0 - pf),
        mutual_info_initial=mi0,
        mutual_info_final=miT,
        mutual_info_retention=mi_ret,
        entropy_final=entropy_vn(rhoT),
        H_stationarity=H_stat,
        B_pop_pair=B_pop,
        B_corr_stat=B_corr_stat,
        S_ent=S_ent,
        V_corr=V_corr,
        R_resp_pop=R_resp_pop,
        R_rob_pop=R_rob_pop,
        F_scalar_pop=F_pop,
        R_resp_stat=R_resp_stat,
        R_rob_stat=R_rob_stat,
        F_scalar_stat=F_stat,
        stationarity_gain=B_corr_stat - B_pop,
        **feats,
    )


def safe_spearman(y: np.ndarray, pred: np.ndarray) -> float:
    rho = spearmanr(y, pred).correlation
    if rho is None or np.isnan(rho):
        return 0.0
    return float(rho)


def evaluate_feature_set(df: pd.DataFrame, name: str, features: list[str], y_col: str, groups: np.ndarray) -> dict[str, float]:
    X = df[features].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    rf = RandomForestRegressor(n_estimators=120, min_samples_leaf=4, random_state=SEED, n_jobs=1)
    pred = cross_val_predict(rf, X, y, cv=kf, n_jobs=1)
    logo = LeaveOneGroupOut()
    pred_logo = cross_val_predict(
        RandomForestRegressor(n_estimators=90, min_samples_leaf=4, random_state=SEED + 1, n_jobs=1),
        X,
        y,
        cv=logo.split(X, y, groups=groups),
        n_jobs=1,
    )
    ridge = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
    pred_ridge = cross_val_predict(ridge, X, y, cv=kf)
    return {
        "feature_set": name,
        "n_features": len(features),
        "rf_r2": float(r2_score(y, pred)),
        "rf_rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "rf_spearman": safe_spearman(y, pred),
        "ridge_r2": float(r2_score(y, pred_ridge)),
        "ridge_spearman": safe_spearman(y, pred_ridge),
        "leave_channel_r2": float(r2_score(y, pred_logo)),
        "leave_channel_spearman": safe_spearman(y, pred_logo),
    }


def scalar_correlations(df: pd.DataFrame, y_col: str) -> pd.DataFrame:
    rows = []
    for col in [
        "H_stationarity",
        "B_pop_pair",
        "B_corr_stat",
        "S_ent",
        "V_corr",
        "R_rob_pop",
        "R_rob_stat",
        "F_scalar_pop",
        "F_scalar_stat",
        "initial_concurrence",
        "dfs_collective_overlap",
        "dfs_parity_overlap",
        "concurrence_slope0",
        "purity_decay0",
    ]:
        sign = -1.0 if col.startswith("F_") or col.endswith("slope0") or col.endswith("decay0") else 1.0
        rows.append({"score": col, "spearman": safe_spearman(df[y_col].to_numpy(), sign * df[col].to_numpy())})
    return pd.DataFrame(rows).sort_values("spearman", ascending=False)


def grouped_importance(df: pd.DataFrame, features: list[str], y_col: str) -> pd.DataFrame:
    X = df[features].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    model = RandomForestRegressor(n_estimators=140, min_samples_leaf=4, random_state=SEED, n_jobs=1)
    model.fit(X, y)
    imp = permutation_importance(model, X, y, n_repeats=6, random_state=SEED, n_jobs=1)
    rows = []
    for feature, mean, std in zip(features, imp.importances_mean, imp.importances_std):
        rows.append({"feature": feature, "importance_mean": float(mean), "importance_std": float(std)})
    return pd.DataFrame(rows).sort_values("importance_mean", ascending=False)


def make_plots(df: pd.DataFrame, comp: pd.DataFrame, corr: pd.DataFrame, imp: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    colors = {
        "rates_only": "#607d8b",
        "einselection_baseline": "#1976d2",
        "short_time_generator_baseline": "#5e35b1",
        "scalar_population": "#ef6c00",
        "scalar_correlated_stationarity": "#00897b",
        "defect_field_stationary": "#2e7d32",
        "shuffled_defect_null": "#9e9e9e",
    }

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    c = comp.sort_values("rf_r2", ascending=False)
    ax.bar(c["feature_set"], c["rf_r2"], color=[colors.get(x, "#777777") for x in c["feature_set"]])
    ax.axhline(0.0, color="#333333", lw=1)
    ax.set_ylabel("5-fold RF R2")
    ax.set_title("Two-qubit entanglement-robustness prediction")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig1_model_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.scatter(df["B_pop_pair"], df["entanglement_robustness"], s=11, alpha=0.35, label="B_pop", color="#ef6c00")
    ax.scatter(df["B_corr_stat"], df["entanglement_robustness"], s=11, alpha=0.35, label="B_corr_stat", color="#00897b")
    ax.set_xlabel("balance support")
    ax.set_ylabel("entanglement robustness target")
    ax.set_title("Population balance versus correlated stationarity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig2_balance_scatter.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ch = df.groupby("channel_family", as_index=False).agg(
        final_concurrence=("final_concurrence", "mean"),
        entanglement_robustness=("entanglement_robustness", "mean"),
        B_corr_stat=("B_corr_stat", "mean"),
    )
    x = np.arange(len(ch))
    width = 0.27
    ax.bar(x - width, ch["entanglement_robustness"], width, label="target", color="#2e7d32")
    ax.bar(x, ch["final_concurrence"], width, label="final concurrence", color="#1976d2")
    ax.bar(x + width, ch["B_corr_stat"], width, label="B_corr_stat", color="#00897b")
    ax.set_xticks(x)
    ax.set_xticklabels(ch["channel_family"], rotation=25, ha="right")
    ax.set_ylabel("mean value")
    ax.set_title("Channel-family structure of entanglement survival")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig3_channel_family_summary.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    top = imp.head(14).iloc[::-1]
    ax.barh(top["feature"], top["importance_mean"], xerr=top["importance_std"], color="#2e7d32", alpha=0.88)
    ax.set_xlabel("permutation importance")
    ax.set_title("Defect-field feature importance")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig4_feature_importance.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    c2 = corr.head(12).iloc[::-1]
    ax.barh(c2["score"], c2["spearman"], color="#00897b")
    ax.set_xlabel("Spearman correlation with target")
    ax.set_title("Scalar support correlations")
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig5_scalar_correlations.png", dpi=180)
    plt.close(fig)


def write_readme(comp: pd.DataFrame, corr: pd.DataFrame) -> None:
    best = comp.sort_values("rf_r2", ascending=False).iloc[0]
    stat = comp[comp["feature_set"] == "scalar_correlated_stationarity"].iloc[0]
    pop = comp[comp["feature_set"] == "scalar_population"].iloc[0]
    text = f"""# Paper 64: Two-Qubit Pointer Selection

This experiment extends the one-qubit pointer benchmarks of Papers 51 and 58
to two qubits. It asks which entangled states survive open-system Lindblad
evolution and whether correlated stationarity is a better balance coordinate
than raw local population balance.

Headline result:

- best feature set: `{best['feature_set']}` with RF R2 `{best['rf_r2']:.4f}`
- scalar correlated stationarity RF R2: `{stat['rf_r2']:.4f}`
- scalar population RF R2: `{pop['rf_r2']:.4f}`
- stationarity gain over population: `{stat['rf_r2'] - pop['rf_r2']:+.4f}`

Outputs:

- `paper64_two_qubit_instances.csv`
- `paper64_model_comparison.csv`
- `paper64_scalar_correlations.csv`
- `paper64_feature_importance.csv`
- `paper64_summary.json`
- `fig1_model_comparison.png`
- `fig2_balance_scatter.png`
- `fig3_channel_family_summary.png`
- `fig4_feature_importance.png`
- `fig5_scalar_correlations.png`

No external quantum dataset is redistributed. All CSV/JSON/PNG files are
derived reproducibility artifacts generated locally by the script.
"""
    (OUTDIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    rng = np.random.default_rng(SEED)
    rows = []
    for idx in range(N_INSTANCES):
        rows.append(asdict(make_record(idx, rng)))
    df = pd.DataFrame(rows)
    y_col = "entanglement_robustness"
    groups = df["channel_family"].to_numpy()

    rates = ["gamma_lz", "gamma_cz", "gamma_zz", "gamma_amp", "gamma_dep", "gamma_heat", "J", "T", "hx1", "hx2", "hz1", "hz2"]
    einselection = [
        "initial_concurrence",
        "bell_phi_overlap",
        "bell_psi_overlap",
        "dfs_collective_overlap",
        "dfs_parity_overlap",
        "local_dephasing_pressure",
        "collective_pressure",
        "pair_dephasing_pressure",
        "damping_pressure",
        "depolarizing_pressure",
        "hamiltonian_coupling_pressure",
    ]
    short_time = ["initial_concurrence", "lindblad_norm0", "concurrence_slope0", "purity_decay0", "local_drift0", "corr_drift0"]
    scalar_pop = ["H_stationarity", "B_pop_pair", "S_ent", "V_corr", "R_rob_pop", "F_scalar_pop"]
    scalar_stat = ["H_stationarity", "B_corr_stat", "S_ent", "V_corr", "R_rob_stat", "F_scalar_stat"]
    defect_field = sorted(
        set(
            rates
            + einselection
            + short_time
            + scalar_stat
            + [
                "coherence_norm0",
                "symmetry_expectation",
                "mutual_info_initial",
                "stationarity_gain",
                "B_pop_pair",
                "R_rob_pop",
                "F_scalar_pop",
            ]
        )
    )
    shuffled = defect_field.copy()
    df_shuf = df.copy()
    for col in shuffled:
        df_shuf[col] = rng.permutation(df_shuf[col].to_numpy())

    feature_sets = [
        ("rates_only", rates, df),
        ("einselection_baseline", einselection, df),
        ("short_time_generator_baseline", short_time, df),
        ("scalar_population", scalar_pop, df),
        ("scalar_correlated_stationarity", scalar_stat, df),
        ("defect_field_stationary", defect_field, df),
        ("shuffled_defect_null", shuffled, df_shuf),
    ]

    comp_rows = []
    for name, features, frame in feature_sets:
        comp_rows.append(evaluate_feature_set(frame, name, features, y_col, groups))
    comp = pd.DataFrame(comp_rows).sort_values("rf_r2", ascending=False)
    corr = scalar_correlations(df, y_col)
    imp = grouped_importance(df, defect_field, y_col)

    df.to_csv(OUTDIR / "paper64_two_qubit_instances.csv", index=False)
    comp.to_csv(OUTDIR / "paper64_model_comparison.csv", index=False)
    corr.to_csv(OUTDIR / "paper64_scalar_correlations.csv", index=False)
    imp.to_csv(OUTDIR / "paper64_feature_importance.csv", index=False)
    make_plots(df, comp, corr, imp)
    write_readme(comp, corr)

    summary = {
        "experiment": "Paper 64 Two-Qubit Pointer Selection",
        "seed": SEED,
        "n_instances": N_INSTANCES,
        "target": y_col,
        "model_comparison": comp.to_dict(orient="records"),
        "top_scalar_correlations": corr.head(12).to_dict(orient="records"),
        "top_feature_importance": imp.head(14).to_dict(orient="records"),
        "channel_family_counts": df["channel_family"].value_counts().to_dict(),
        "state_family_counts": df["state_family"].value_counts().to_dict(),
    }
    (OUTDIR / "paper64_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Paper 64 two-qubit pointer benchmark complete.")
    print(comp.to_string(index=False))
    print()
    print(corr.head(12).to_string(index=False))
    print(f"\nOutputs written to: {OUTDIR.resolve()}")


if __name__ == "__main__":
    main()

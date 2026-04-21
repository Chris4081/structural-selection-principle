#!/usr/bin/env python3
"""Small ensemble scan for the reduced type IIB / KKLT structural master functional."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import root


@dataclass(frozen=True)
class DynamicParams:
    W0: float
    A: float
    a: float
    D: float
    n: float


@dataclass(frozen=True)
class StructuralParams:
    alpha_h: float = 1.0
    alpha_b: float = 1.0
    alpha_s: float = 1.0
    alpha_r: float = 1.0
    beta_r: float = 1.0
    alpha_sigma: float = 1.0
    beta_sigma: float = 1.0
    gamma_sigma: float = 1.0
    epsilon_v: float = 1.0e-6
    tau_control: float = 10.0
    xi_hat: float = 1.0
    tadpole_ref: float = 16.0
    m_star_sq: float = 1.0e-14
    weights: tuple[float, float, float, float, float, float] = (1 / 6,) * 6


def parse_csv_floats(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def normalize_theta(theta: float, a_param: float) -> float:
    period = 2.0 * math.pi / a_param
    theta = math.fmod(theta, period)
    if theta < 0.0:
        theta += period
    if theta > period / 2.0:
        theta -= period
    return theta


def kklt_potential(tau: float, theta: float, params: DynamicParams) -> float:
    if tau <= 0.0:
        return float("inf")
    t_field = tau + 1j * theta
    two_tau = 2.0 * tau
    k_val = -3.0 * math.log(two_tau)
    w_val = params.W0 + params.A * np.exp(-params.a * t_field)
    dw_val = -params.a * params.A * np.exp(-params.a * t_field)
    dk_val = -3.0 / two_tau
    dtw_val = dw_val + dk_val * w_val
    kinv_val = two_tau**2 / 3.0
    vf_val = math.exp(k_val) * (kinv_val * abs(dtw_val) ** 2 - 3.0 * abs(w_val) ** 2)
    vup_val = params.D / (two_tau**params.n)
    return float(np.real(vf_val) + vup_val)


def central_gradient(tau: float, theta: float, params: DynamicParams) -> np.ndarray:
    h_tau = 1.0e-5 * max(1.0, abs(tau))
    h_theta = 1.0e-5 * max(1.0, abs(theta), 1.0)
    f = lambda x_tau, x_theta: kklt_potential(x_tau, x_theta, params)
    d_tau = (f(tau + h_tau, theta) - f(tau - h_tau, theta)) / (2.0 * h_tau)
    d_theta = (f(tau, theta + h_theta) - f(tau, theta - h_theta)) / (2.0 * h_theta)
    return np.array([d_tau, d_theta], dtype=float)


def central_hessian(tau: float, theta: float, params: DynamicParams) -> np.ndarray:
    h_tau = 1.0e-4 * max(1.0, abs(tau))
    h_theta = 1.0e-4 * max(1.0, abs(theta), 1.0)
    f = lambda x_tau, x_theta: kklt_potential(x_tau, x_theta, params)
    f00 = f(tau, theta)
    dtt = (f(tau + h_tau, theta) - 2.0 * f00 + f(tau - h_tau, theta)) / (h_tau**2)
    dyy = (f(tau, theta + h_theta) - 2.0 * f00 + f(tau, theta - h_theta)) / (h_theta**2)
    dty = (
        f(tau + h_tau, theta + h_theta)
        - f(tau + h_tau, theta - h_theta)
        - f(tau - h_tau, theta + h_theta)
        + f(tau - h_tau, theta - h_theta)
    ) / (4.0 * h_tau * h_theta)
    return np.array([[dtt, dty], [dty, dyy]], dtype=float)


def canonical_metric_factor(tau: float) -> float:
    return 3.0 / (4.0 * tau * tau)


def canonical_hessian(tau: float, theta: float, params: DynamicParams) -> np.ndarray:
    return central_hessian(tau, theta, params) / canonical_metric_factor(tau)


def solve_stationary_points(
    params: DynamicParams,
    tau_seeds: Iterable[float],
    theta_seeds: Iterable[float],
    tau_min: float,
    tau_max: float,
    grad_tol: float,
) -> list[tuple[float, float]]:
    solutions: list[tuple[float, float]] = []

    def equations(x: np.ndarray) -> np.ndarray:
        tau, theta = float(x[0]), float(x[1])
        if tau <= tau_min / 10.0:
            return np.array([1.0e3, 1.0e3], dtype=float)
        return central_gradient(tau, theta, params)

    for tau0, theta0 in product(tau_seeds, theta_seeds):
        result = root(equations, np.array([tau0, theta0], dtype=float), method="hybr")
        if not result.success:
            continue
        tau, theta = float(result.x[0]), normalize_theta(float(result.x[1]), params.a)
        if not (tau_min <= tau <= tau_max):
            continue
        if np.linalg.norm(central_gradient(tau, theta, params)) > grad_tol:
            continue
        duplicate = False
        for tau_old, theta_old in solutions:
            if abs(tau - tau_old) < 1.0e-4 and abs(theta - theta_old) < 1.0e-4:
                duplicate = True
                break
        if not duplicate:
            solutions.append((tau, theta))
    return sorted(solutions)


def classify_stationary_point(
    energy: float,
    lambda_min: float,
    stability_tol: float = 1.0e-16,
    energy_tol: float = 1.0e-16,
) -> str:
    if lambda_min > stability_tol:
        if energy > energy_tol:
            return "dS minimum"
        if abs(energy) <= energy_tol:
            return "near-Minkowski minimum"
        return "AdS minimum"
    return "unstable stationary point"


def estimate_decay_rate(
    tau: float,
    theta: float,
    energy: float,
    params: DynamicParams,
    tau_barrier_max: float,
    barrier_samples: int,
) -> tuple[float, float]:
    tau_stop = max(tau_barrier_max, tau * 1.05)
    taus = np.linspace(tau, tau_stop, barrier_samples)
    values = np.array([kklt_potential(x_tau, theta, params) for x_tau in taus], dtype=float)
    barrier = max(0.0, float(np.max(values) - energy))
    scale = max(abs(energy), np.max(np.abs(values)), abs(params.W0) ** 2, 1.0e-18)
    if barrier <= 0.0:
        return 1.0, 0.0
    exponent = min(barrier / scale, 700.0)
    return float(math.exp(-exponent)), barrier


def candidate_record(
    params: DynamicParams,
    structural: StructuralParams,
    delta_q_d3: float,
    tau: float,
    theta: float,
    tau_barrier_max: float,
    barrier_samples: int,
) -> dict:
    energy = kklt_potential(tau, theta, params)
    grad = central_gradient(tau, theta, params)
    hessian = canonical_hessian(tau, theta, params)
    eigvals = np.linalg.eigvalsh(hessian)
    lambda_min = float(np.min(eigvals))
    vacuum_type = classify_stationary_point(energy, lambda_min)

    v_scale = max(abs(params.W0) ** 2, abs(energy), 1.0e-18)
    d_tau_hat = float(grad[0] / v_scale)
    d_theta_hat = float(grad[1] / v_scale)

    residual_h = structural.alpha_h * (
        d_tau_hat**2 + (d_theta_hat**2) / max(tau * tau, 1.0e-12)
    )

    residual_b = structural.alpha_b * (delta_q_d3 / structural.tadpole_ref) ** 2

    sigma_vac = (
        d_tau_hat**2 + (d_theta_hat**2) / max(tau * tau, 1.0e-12)
    ) / (
        1.0 + d_tau_hat**2 + (d_theta_hat**2) / max(tau * tau, 1.0e-12)
    )
    residual_s = structural.alpha_s * sigma_vac

    offdiag = float(abs(hessian[0, 1]))
    diag_a = float(abs(hessian[0, 0]) + structural.m_star_sq)
    diag_b = float(abs(hessian[1, 1]) + structural.m_star_sq)
    interaction_v = offdiag / math.sqrt(diag_a * diag_b)
    score_v = interaction_v / (1.0 + interaction_v)
    v_epsilon = structural.epsilon_v + (1.0 - structural.epsilon_v) * score_v

    decay_rate, barrier_height = estimate_decay_rate(
        tau=tau,
        theta=theta,
        energy=energy,
        params=params,
        tau_barrier_max=tau_barrier_max,
        barrier_samples=barrier_samples,
    )
    mass_ref_sq = max(float(np.max(np.abs(eigvals))), 1.0e-18)
    residual_r = structural.alpha_r * max(0.0, -lambda_min / mass_ref_sq) + structural.beta_r * decay_rate

    residual_sigma = (
        structural.alpha_sigma * (structural.tau_control / tau) ** 2
        + structural.beta_sigma * math.exp(-params.a * tau)
        + structural.gamma_sigma * (structural.xi_hat / ((2.0 * tau) ** 1.5)) ** 2
    )

    weights = structural.weights
    f_master = (
        weights[0] * residual_h
        + weights[1] * residual_b
        + weights[2] * residual_s
        + weights[3] * (-math.log(v_epsilon))
        + weights[4] * residual_r
        + weights[5] * residual_sigma
    )

    scores = {
        "H": math.exp(-residual_h),
        "B": math.exp(-residual_b),
        "S": math.exp(-residual_s),
        "V": score_v,
        "R": math.exp(-residual_r),
        "Sigma": math.exp(-residual_sigma),
    }

    return {
        "params": asdict(params),
        "delta_q_d3": delta_q_d3,
        "tau": tau,
        "theta": theta,
        "energy": energy,
        "vacuum_type": vacuum_type,
        "lambda_min": lambda_min,
        "hessian_eigenvalues": eigvals.tolist(),
        "barrier_height_tau": barrier_height,
        "decay_rate_heuristic": decay_rate,
        "residuals": {
            "R_H": residual_h,
            "R_B": residual_b,
            "R_S": residual_s,
            "R_R": residual_r,
            "R_Sigma": residual_sigma,
            "minus_log_Veps": -math.log(v_epsilon),
        },
        "scores": scores,
        "connectivity_aux": {
            "I_V": interaction_v,
            "V_epsilon": v_epsilon,
        },
        "f_master": f_master,
    }


def summarize(candidates: list[dict]) -> dict:
    by_master = sorted(candidates, key=lambda item: item["f_master"])
    by_energy = sorted(candidates, key=lambda item: item["energy"])
    class_counts: dict[str, int] = {}
    for item in candidates:
        class_counts[item["vacuum_type"]] = class_counts.get(item["vacuum_type"], 0) + 1
    return {
        "num_candidates": len(candidates),
        "class_counts": class_counts,
        "top_by_master": by_master[:5],
        "top_by_energy": by_energy[:5],
    }


def write_json(path: Path, payload: dict, pretty: bool) -> None:
    path.write_text(
        json.dumps(payload, indent=2 if pretty else None, sort_keys=pretty),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan a reduced KKLT ensemble and evaluate the structural master functional."
    )
    parser.add_argument("--w0-values", default="-5e-5,-1e-4,-2e-4,-3e-4")
    parser.add_argument("--d-values", default="0,1e-9")
    parser.add_argument("--delta-q-values", default="0,8")
    parser.add_argument("--A", type=float, default=1.0)
    parser.add_argument("--a", type=float, default=0.1)
    parser.add_argument("--n", type=float, default=2.0)
    parser.add_argument("--tau-seed-min", type=float, default=20.0)
    parser.add_argument("--tau-seed-max", type=float, default=160.0)
    parser.add_argument("--num-tau-seeds", type=int, default=18)
    parser.add_argument("--tau-max-stationary", type=float, default=160.0)
    parser.add_argument("--tau-min-stationary", type=float, default=5.0)
    parser.add_argument("--tau-barrier-max", type=float, default=220.0)
    parser.add_argument("--barrier-samples", type=int, default=300)
    parser.add_argument("--grad-tol", type=float, default=1.0e-8)
    parser.add_argument("--json-output")
    parser.add_argument("--pretty-json", action="store_true")
    args = parser.parse_args()

    structural = StructuralParams()
    dynamic_backgrounds = [
        DynamicParams(W0=w0, A=args.A, a=args.a, D=d_value, n=args.n)
        for w0, d_value in product(parse_csv_floats(args.w0_values), parse_csv_floats(args.d_values))
    ]
    delta_q_values = parse_csv_floats(args.delta_q_values)

    tau_seeds = np.linspace(args.tau_seed_min, args.tau_seed_max, args.num_tau_seeds)
    theta_seeds = [0.0, math.pi / args.a]

    candidates: list[dict] = []
    dynamic_solutions: list[dict] = []
    for params in dynamic_backgrounds:
        points = solve_stationary_points(
            params=params,
            tau_seeds=tau_seeds,
            theta_seeds=theta_seeds,
            tau_min=args.tau_min_stationary,
            tau_max=args.tau_max_stationary,
            grad_tol=args.grad_tol,
        )
        dynamic_solutions.append(
            {
                "params": asdict(params),
                "stationary_points": [{"tau": tau, "theta": theta} for tau, theta in points],
            }
        )
        for tau, theta in points:
            for delta_q in delta_q_values:
                candidates.append(
                    candidate_record(
                        params=params,
                        structural=structural,
                        delta_q_d3=delta_q,
                        tau=tau,
                        theta=theta,
                        tau_barrier_max=args.tau_barrier_max,
                        barrier_samples=args.barrier_samples,
                    )
                )

    summary = summarize(candidates)
    payload = {
        "scan_config": {
            "dynamic_backgrounds": [asdict(item) for item in dynamic_backgrounds],
            "delta_q_values": delta_q_values,
            "tau_seed_min": args.tau_seed_min,
            "tau_seed_max": args.tau_seed_max,
            "num_tau_seeds": args.num_tau_seeds,
            "tau_max_stationary": args.tau_max_stationary,
            "tau_barrier_max": args.tau_barrier_max,
            "grad_tol": args.grad_tol,
        },
        "structural_params": asdict(structural),
        "dynamic_solutions": dynamic_solutions,
        "candidates": candidates,
        "summary": summary,
    }

    print(f"Scanned {len(dynamic_backgrounds)} dynamic backgrounds.")
    print(f"Found {sum(len(item['stationary_points']) for item in dynamic_solutions)} stationary points.")
    print(f"Constructed {len(candidates)} structural candidates.")
    print("Top candidates by F_master:")
    for rank, item in enumerate(summary["top_by_master"], start=1):
        params = item["params"]
        print(
            f"  {rank}. F={item['f_master']:.6f} | {item['vacuum_type']} | "
            f"W0={params['W0']:.3e}, D={params['D']:.3e}, "
            f"tau={item['tau']:.3f}, theta={item['theta']:.3f}, "
            f"deltaQ={item['delta_q_d3']:.1f}, V={item['energy']:.3e}"
        )

    if args.json_output:
        write_json(Path(args.json_output), payload, args.pretty_json)
        print(f"Wrote JSON results to {args.json_output}")


if __name__ == "__main__":
    main()

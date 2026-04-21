#!/usr/bin/env python3
"""Numerical benchmark for the static phi^4 structural-selection protocol.

This script implements the lattice diagnostics and benchmark tests introduced in
``structural_selection_phi4.tex``:

1. Vacuum discrimination in the trivial sector.
2. Kink discrimination in the topological sector.
3. Structural ranking versus static-energy ranking.

The implementation stays close to the paper definitions while adding one useful
numerical refinement for the kink sector: the translation-like mode is detected
by overlap with the discrete derivative of the profile, rather than by blindly
dropping the second-smallest Hessian eigenvalue for every q=1 configuration.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scipy.linalg import eigh_tridiagonal
except ImportError:  # pragma: no cover - fallback path
    eigh_tridiagonal = None


@dataclass
class ProtocolConfig:
    v: float = 1.0
    lambda_phi: float = 0.5
    L: float = 20.0
    n_x: int = 801
    ell_s: float = 1.0
    ell_v: float = 1.0
    gamma_r: float = 4.0
    alpha_h: float = 4.0
    alpha_b: float = 1.0
    alpha_s: float = 0.0
    alpha_v: float = 1.0
    alpha_r: float = 3.0
    alpha_q: float = 10.0
    sigma: float = 5.0
    relaxation_dt: float = 8.0e-4
    relaxation_tol: float = 1.0e-10
    relaxation_max_steps: int = 200000
    moderate_eta_max: float = 0.30
    kink_success_threshold: float = 0.90


class Phi4StaticProtocol:
    def __init__(self, config: ProtocolConfig) -> None:
        if config.n_x < 5 or config.n_x % 2 == 0:
            raise ValueError("n_x must be an odd integer >= 5")

        self.cfg = config
        self.a = 2.0 * config.L / (config.n_x - 1)
        self.x = np.linspace(-config.L, config.L, config.n_x)
        self.mass_sq = 2.0 * config.lambda_phi * config.v**2
        self.n_h = 2.0 * config.L
        self.n_b = 2.0 * config.L
        self.n_s = 2.0 * config.L
        self.n_v = (2.0 * config.L) ** 2

        delta = self.x[:, None] - self.x[None, :]
        self.kernel = (
            np.exp(-(delta**2) / (2.0 * config.ell_v**2))
            / (math.sqrt(2.0 * math.pi) * config.ell_v)
        )
        self.kernel_row_sums = self.kernel.sum(axis=1)

    @property
    def interior(self) -> slice:
        return slice(1, -1)

    def gradient(self, phi: np.ndarray) -> np.ndarray:
        return (phi[2:] - phi[:-2]) / (2.0 * self.a)

    def laplacian(self, phi: np.ndarray) -> np.ndarray:
        return (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / (self.a**2)

    def residual(self, phi: np.ndarray) -> np.ndarray:
        interior_phi = phi[1:-1]
        return -self.laplacian(phi) + self.cfg.lambda_phi * interior_phi * (
            interior_phi**2 - self.cfg.v**2
        )

    def potential(self, phi: np.ndarray) -> np.ndarray:
        return 0.25 * self.cfg.lambda_phi * (phi**2 - self.cfg.v**2) ** 2

    def topological_charge(self, phi: np.ndarray) -> float:
        return (phi[-1] - phi[0]) / (2.0 * self.cfg.v)

    def harmony(self, phi: np.ndarray) -> tuple[float, float]:
        q_h = (self.a / self.n_h) * float(np.sum(self.residual(phi) ** 2))
        return math.exp(-q_h), q_h

    def balance(self, phi: np.ndarray) -> tuple[float, float]:
        b = self.residual(phi) * self.gradient(phi)
        q_b = (self.a / self.n_b) * float(np.sum(b**2))
        return math.exp(-q_b), q_b

    def activity(self, phi: np.ndarray) -> tuple[float, float]:
        q_s = (
            self.a
            * self.cfg.ell_s**2
            / self.n_s
            * float(np.sum(self.gradient(phi) ** 2))
        )
        s = q_s / (1.0 + q_s)
        return s, q_s

    def connectivity(self, phi: np.ndarray) -> tuple[float, float]:
        quad = float(phi @ (self.kernel @ phi))
        weighted_norm = float((phi**2) @ self.kernel_row_sums)
        q_v = (2.0 * self.a**2 / self.n_v) * (weighted_norm - quad)
        q_v = max(q_v, 0.0)
        return math.exp(-q_v), q_v

    def _hessian_tridiagonal(self, phi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        interior_phi = phi[1:-1]
        diag = 2.0 / self.a**2 + self.cfg.lambda_phi * (
            3.0 * interior_phi**2 - self.cfg.v**2
        )
        offdiag = np.full(interior_phi.size - 1, -1.0 / self.a**2)
        return diag, offdiag

    def _lowest_eigensystem(
        self, diag: np.ndarray, offdiag: np.ndarray, count: int
    ) -> tuple[np.ndarray, np.ndarray | None]:
        count = max(1, min(count, diag.size))
        if eigh_tridiagonal is not None:
            eigvals, eigvecs = eigh_tridiagonal(
                diag,
                offdiag,
                select="i",
                select_range=(0, count - 1),
                eigvals_only=False,
            )
            return eigvals, eigvecs

        dense = np.diag(diag) + np.diag(offdiag, k=1) + np.diag(offdiag, k=-1)
        eigvals, eigvecs = np.linalg.eigh(dense)
        return eigvals[:count], eigvecs[:, :count]

    def stability_eigenvalue(self, phi: np.ndarray, q: int) -> float:
        diag, offdiag = self._hessian_tridiagonal(phi)

        if q == 0:
            eigvals, _ = self._lowest_eigensystem(diag, offdiag, count=1)
            return float(eigvals[0])

        eigvals, eigvecs = self._lowest_eigensystem(diag, offdiag, count=6)
        if eigvecs is None:
            return float(eigvals[min(1, len(eigvals) - 1)])

        translation = self.gradient(phi)
        norm = float(np.linalg.norm(translation))
        if norm < 1.0e-14:
            return float(eigvals[min(1, len(eigvals) - 1)])

        translation /= norm
        overlaps = np.abs(eigvecs.T @ translation)
        translation_index = int(np.argmax(overlaps))

        candidates = [
            float(eigvals[i]) for i in range(len(eigvals)) if i != translation_index
        ]
        if not candidates:
            return float(eigvals[-1])
        return min(candidates)

    def robustness(self, phi: np.ndarray, q: int) -> tuple[float, float]:
        lambda_star = self.stability_eigenvalue(phi, q)
        instability = max(0.0, -lambda_star / self.mass_sq)
        r = math.exp(-self.cfg.gamma_r * instability)
        return r, lambda_star

    def static_energy(self, phi: np.ndarray) -> float:
        grad = self.gradient(phi)
        pot = self.potential(phi[1:-1])
        return self.a * float(np.sum(0.5 * grad**2 + pot))

    def evaluate(self, phi: np.ndarray, q: int) -> dict[str, float]:
        h, q_h = self.harmony(phi)
        b, q_b = self.balance(phi)
        s, q_s = self.activity(phi)
        v_val, q_v = self.connectivity(phi)
        r_val, lambda_star = self.robustness(phi, q=q)
        q_top = float(self.topological_charge(phi))

        f_struct = float(
            self.cfg.alpha_h * (1.0 - h)
            + self.cfg.alpha_b * (1.0 - b)
            + self.cfg.alpha_s * s
            + self.cfg.alpha_v * (1.0 - v_val)
            + self.cfg.alpha_r * (1.0 - r_val)
            + self.cfg.alpha_q * (q_top - q) ** 2
        )

        return {
            "H": float(h),
            "Q_H": float(q_h),
            "B": float(b),
            "Q_B": float(q_b),
            "S": float(s),
            "Q_S": float(q_s),
            "V": float(v_val),
            "Q_V": float(q_v),
            "R": float(r_val),
            "lambda_star": float(lambda_star),
            "Q_top": q_top,
            "F_struct": f_struct,
            "E_stat": float(self.static_energy(phi)),
        }

    def vacuum_state(self, sign: float) -> np.ndarray:
        return np.full(self.cfg.n_x, sign * self.cfg.v, dtype=float)

    def saddle_state(self) -> np.ndarray:
        return np.zeros(self.cfg.n_x, dtype=float)

    def kink_initial_guess(self) -> np.ndarray:
        phi = self.cfg.v * np.tanh(0.5 * math.sqrt(self.mass_sq) * self.x)
        phi[0] = -self.cfg.v
        phi[-1] = self.cfg.v
        return phi

    def relax_kink(self) -> tuple[np.ndarray, dict[str, float]]:
        phi = self.kink_initial_guess()
        meta: dict[str, float] = {"steps": 0.0, "max_update": 0.0, "max_residual": 0.0}

        for step in range(1, self.cfg.relaxation_max_steps + 1):
            residual = self.residual(phi)
            update = -self.cfg.relaxation_dt * residual
            phi[1:-1] += update
            phi[0] = -self.cfg.v
            phi[-1] = self.cfg.v

            max_update = float(np.max(np.abs(update)))
            max_residual = float(np.max(np.abs(residual)))
            meta = {
                "steps": float(step),
                "max_update": max_update,
                "max_residual": max_residual,
            }
            if max_update < self.cfg.relaxation_tol and max_residual < self.cfg.relaxation_tol:
                break

        return phi, meta

    def perturbed_vacuum(self, sign: float, eta: float) -> np.ndarray:
        phi = sign * self.cfg.v + eta * np.cos(math.pi * self.x / (2.0 * self.cfg.L))
        phi[0] = sign * self.cfg.v
        phi[-1] = sign * self.cfg.v
        return phi

    def distorted_kink(self, kink: np.ndarray, epsilon: float, mode: int) -> np.ndarray:
        envelope = np.exp(-(self.x**2) / (self.cfg.sigma**2))
        sinusoid = np.sin(mode * math.pi * (self.x + self.cfg.L) / (2.0 * self.cfg.L))
        phi = kink + epsilon * envelope * sinusoid
        phi[0] = -self.cfg.v
        phi[-1] = self.cfg.v
        return phi

    def kink_distance(self, phi: np.ndarray, kink: np.ndarray) -> float:
        return math.sqrt(self.a / (2.0 * self.cfg.L) * float(np.sum((phi - kink) ** 2)))


def protocol_a(
    protocol: Phi4StaticProtocol,
) -> dict[str, Any]:
    plus = protocol.vacuum_state(+1.0)
    minus = protocol.vacuum_state(-1.0)
    saddle = protocol.saddle_state()

    plus_eval = protocol.evaluate(plus, q=0)
    minus_eval = protocol.evaluate(minus, q=0)
    saddle_eval = protocol.evaluate(saddle, q=0)

    perturbation_rows = []
    eta_values = [round(0.05 * i, 2) for i in range(1, 11)]
    moderate_passes = []
    for eta in eta_values:
        state = protocol.perturbed_vacuum(+1.0, eta)
        evaluation = protocol.evaluate(state, q=0)
        ordering_pass = (
            plus_eval["F_struct"] < evaluation["F_struct"] < saddle_eval["F_struct"]
        )
        perturbation_rows.append(
            {
                "eta": eta,
                "F_struct": evaluation["F_struct"],
                "ordering_pass": ordering_pass,
            }
        )
        if eta <= protocol.cfg.moderate_eta_max:
            moderate_passes.append(ordering_pass)

    delta_vac = saddle_eval["F_struct"] - plus_eval["F_struct"]

    return {
        "exact_states": {
            "vacuum_plus": plus_eval,
            "vacuum_minus": minus_eval,
            "saddle": saddle_eval,
        },
        "delta_vac": delta_vac,
        "vacuum_vs_saddle_pass": plus_eval["F_struct"] < saddle_eval["F_struct"]
        and minus_eval["F_struct"] < saddle_eval["F_struct"],
        "moderate_eta_max": protocol.cfg.moderate_eta_max,
        "moderate_perturbation_pass": all(moderate_passes),
        "perturbations": perturbation_rows,
    }


def protocol_b(
    protocol: Phi4StaticProtocol, kink: np.ndarray
) -> dict[str, Any]:
    kink_eval = protocol.evaluate(kink, q=1)

    distortion_rows = []
    better_count = 0
    total = 0
    for mode in (1, 2, 3):
        for epsilon in [round(0.02 * i, 2) for i in range(1, 16)]:
            state = protocol.distorted_kink(kink, epsilon, mode)
            evaluation = protocol.evaluate(state, q=1)
            better = evaluation["F_struct"] > kink_eval["F_struct"]
            better_count += int(better)
            total += 1
            distortion_rows.append(
                {
                    "mode": mode,
                    "epsilon": epsilon,
                    "F_struct": evaluation["F_struct"],
                    "E_stat": evaluation["E_stat"],
                    "distance_to_kink": protocol.kink_distance(state, kink),
                    "kink_outperforms": better,
                }
            )

    p_k = better_count / total if total else 0.0
    return {
        "kink": kink_eval,
        "distortions": distortion_rows,
        "p_k": p_k,
        "success_threshold": protocol.cfg.kink_success_threshold,
        "pass": p_k >= protocol.cfg.kink_success_threshold,
    }


def protocol_c(
    protocol_b_results: dict[str, Any],
) -> dict[str, Any]:
    distortions = protocol_b_results["distortions"]
    percentiles = [5, 10, 20]

    by_f = sorted(distortions, key=lambda row: row["F_struct"])
    by_e = sorted(distortions, key=lambda row: row["E_stat"])

    results = []
    n = len(distortions)
    for percentile in percentiles:
        count = max(1, math.ceil(percentile / 100.0 * n))
        mean_f = sum(row["distance_to_kink"] for row in by_f[:count]) / count
        mean_e = sum(row["distance_to_kink"] for row in by_e[:count]) / count
        results.append(
            {
                "percentile": percentile,
                "count": count,
                "mean_distance_structural": mean_f,
                "mean_distance_energy": mean_e,
                "structural_better": mean_f < mean_e,
            }
        )

    return {
        "rank_comparison": results,
        "any_structural_gain": any(row["structural_better"] for row in results),
    }


def run_protocol(config: ProtocolConfig) -> dict[str, Any]:
    protocol = Phi4StaticProtocol(config)
    kink, relaxation = protocol.relax_kink()

    a_results = protocol_a(protocol)
    b_results = protocol_b(protocol, kink)
    c_results = protocol_c(b_results)

    return {
        "config": {
            **asdict(config),
            "a": protocol.a,
            "mass_sq": protocol.mass_sq,
        },
        "kink_relaxation": relaxation,
        "protocol_a": a_results,
        "protocol_b": b_results,
        "protocol_c": c_results,
    }


def to_builtin(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: to_builtin(item) for key, item in value.items()}
    if isinstance(value, list):
        return [to_builtin(item) for item in value]
    if isinstance(value, tuple):
        return [to_builtin(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def print_summary(results: dict[str, Any]) -> None:
    a_results = results["protocol_a"]
    b_results = results["protocol_b"]
    c_results = results["protocol_c"]

    vac = a_results["exact_states"]
    print("Protocol A: vacuum discrimination")
    print(
        f"  F(vac+)={vac['vacuum_plus']['F_struct']:.6f}  "
        f"F(vac-)={vac['vacuum_minus']['F_struct']:.6f}  "
        f"F(saddle)={vac['saddle']['F_struct']:.6f}"
    )
    print(
        f"  delta_vac={a_results['delta_vac']:.6f}  "
        f"vacuum_vs_saddle_pass={a_results['vacuum_vs_saddle_pass']}  "
        f"moderate_perturbation_pass={a_results['moderate_perturbation_pass']}"
    )

    print("\nProtocol B: kink discrimination")
    print(
        f"  F(kink)={b_results['kink']['F_struct']:.6f}  "
        f"p_K={b_results['p_k']:.6f}  "
        f"threshold={b_results['success_threshold']:.2f}  "
        f"pass={b_results['pass']}"
    )

    print("\nProtocol C: structural ranking versus energy ranking")
    for row in c_results["rank_comparison"]:
        print(
            f"  p={row['percentile']:>2}%  n={row['count']:>2}  "
            f"d_F={row['mean_distance_structural']:.6f}  "
            f"d_E={row['mean_distance_energy']:.6f}  "
            f"structural_better={row['structural_better']}"
        )
    print(f"  any_structural_gain={c_results['any_structural_gain']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the static phi^4 structural-selection benchmark."
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="Optional path for a JSON dump of all benchmark results.",
    )
    parser.add_argument(
        "--pretty-json",
        action="store_true",
        help="Pretty-print JSON when writing --json-output.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = run_protocol(ProtocolConfig())
    print_summary(results)

    if args.json_output is not None:
        args.json_output.write_text(
            json.dumps(to_builtin(results), indent=2 if args.pretty_json else None),
            encoding="utf-8",
        )
        print(f"\nWrote JSON results to {args.json_output}")


if __name__ == "__main__":
    main()

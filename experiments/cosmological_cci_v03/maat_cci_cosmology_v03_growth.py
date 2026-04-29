#!/usr/bin/env python3
"""Cosmological CCI v0.3: H(z) + fσ8 growth-connectivity benchmark.

This script extends the v0.2 chronometer-only diagnostic with two explicitly
measured sectors:

1. V_corr(z): growth/connectivity support from fσ8 residuals.
2. R_robust(z): joint H(z)+fσ8 consistency support.

It also derives a data-driven transition redshift z_c from the curvature of
log CCI_v03 and fits companion sector weights lambda_a with a small MaxEnt
matching problem. The benchmark is diagnostic, not a precision likelihood.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/codex-mpl-cache")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/codex-mpl-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.optimize import minimize
except Exception:  # pragma: no cover - fallback for minimal environments.
    minimize = None


EPS = 1.0e-9


@dataclass(frozen=True)
class Config:
    H0: float = 67.4
    omega_m: float = 0.315
    omega_l: float = 0.685
    sigma8_0: float = 0.811
    gamma_growth: float = 0.55
    z_min: float = 0.05
    z_max: float = 2.0
    n_grid: int = 450
    kernel_h: float = 0.12
    kernel_growth: float = 0.08
    maxent_best_fraction: float = 0.20
    maxent_prior_strength: float = 0.02
    maxent_balance_strength: float = 0.01


def e_lcdm(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    return np.sqrt(cfg.omega_m * (1.0 + z) ** 3 + cfg.omega_l)


def omega_m_z(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    ez2 = e_lcdm(z, cfg) ** 2
    return cfg.omega_m * (1.0 + z) ** 3 / ez2


def omega_l_z(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    ez2 = e_lcdm(z, cfg) ** 2
    return cfg.omega_l / ez2


def growth_suppression(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    om = omega_m_z(z, cfg)
    ol = omega_l_z(z, cfg)
    denominator = om ** (4.0 / 7.0) - ol + (1.0 + om / 2.0) * (1.0 + ol / 70.0)
    return (5.0 * om / 2.0) / denominator


def growth_factor(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    return growth_suppression(z, cfg) / (growth_suppression(0.0, cfg) * (1.0 + z))


def growth_rate(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    return omega_m_z(z, cfg) ** cfg.gamma_growth


def fsigma8_ref(z: np.ndarray | float, cfg: Config) -> np.ndarray | float:
    return growth_rate(z, cfg) * cfg.sigma8_0 * growth_factor(z, cfg)


def read_csv_columns(path: Path) -> dict[str, np.ndarray | list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result: dict[str, list] = {key: [] for key in rows[0].keys()}
    for row in rows:
        for key, value in row.items():
            try:
                result[key].append(float(value))
            except ValueError:
                result[key].append(value)
    return {
        key: np.array(values, dtype=float) if isinstance(values[0], float) else values
        for key, values in result.items()
    }


def local_chi2(
    z_grid: np.ndarray,
    z_data: np.ndarray,
    residual_pull: np.ndarray,
    width: float,
) -> np.ndarray:
    dz = z_grid[:, None] - z_data[None, :]
    weights = np.exp(-0.5 * (dz / width) ** 2)
    norm = np.sum(weights, axis=1)
    return np.sum(weights * residual_pull[None, :] ** 2, axis=1) / (norm + EPS)


def cci_v02(z: np.ndarray, cfg: Config) -> np.ndarray:
    e = e_lcdm(z, cfg)
    d = growth_factor(z, cfg)
    u = np.abs(e - d) / (e + d + EPS)
    return e * (1.0 + z) * (1.0 + u) / (1.0 + d + EPS)


def cci_v03(z: np.ndarray, cfg: Config, v_corr: np.ndarray, r_robust: np.ndarray) -> np.ndarray:
    e = e_lcdm(z, cfg)
    d = growth_factor(z, cfg)
    u = np.abs(e - d) / (e + d + EPS)
    return e * (1.0 + z) * (1.0 + u) / (1.0 + d + v_corr + r_robust + EPS)


def transition_redshift(z: np.ndarray, cci_norm: np.ndarray) -> tuple[float, np.ndarray]:
    log_curve = np.log(cci_norm + EPS)
    first = np.gradient(log_curve, z)
    curvature = np.gradient(first, z)
    mask = (z > z.min() + 0.12) & (z < z.max() - 0.12)
    idx_local = np.argmax(np.abs(curvature[mask]))
    idx = np.where(mask)[0][idx_local]
    return float(z[idx]), curvature


def fit_maxent_lambdas(defects: np.ndarray, best_fraction: float, prior_strength: float, balance_strength: float) -> dict:
    sectors = ["H", "S", "V", "R"]
    total = np.sum(defects, axis=1)
    n_best = max(3, int(best_fraction * len(total)))
    best_idx = np.argsort(total)[:n_best]
    targets = np.mean(defects[best_idx], axis=0)
    lambda0 = np.ones(defects.shape[1])

    def loss(log_lambda: np.ndarray) -> float:
        lam = np.exp(log_lambda)
        energy = defects @ lam
        weights = np.exp(-(energy - np.min(energy)))
        weights = weights / (np.sum(weights) + EPS)
        expected = weights @ defects
        shares = lam / (np.sum(lam) + EPS)
        mismatch = np.sum((expected - targets) ** 2)
        prior = prior_strength * np.sum((log_lambda - np.log(lambda0)) ** 2)
        balance = balance_strength * np.sum((shares - 1.0 / len(lam)) ** 2)
        return float(mismatch + prior + balance)

    if minimize is not None:
        result = minimize(loss, np.zeros(defects.shape[1]), method="Nelder-Mead", options={"maxiter": 20000})
        log_lam = result.x
        success = bool(result.success)
        fit_loss = float(result.fun)
    else:
        candidates = [np.zeros(defects.shape[1])]
        rng = np.random.default_rng(42)
        candidates.extend(rng.normal(0.0, 0.6, size=(5000, defects.shape[1])))
        losses = np.array([loss(c) for c in candidates])
        log_lam = candidates[int(np.argmin(losses))]
        success = False
        fit_loss = float(np.min(losses))

    lam = np.exp(log_lam)
    shares = lam / (np.sum(lam) + EPS)
    return {
        "sectors": sectors,
        "lambda": {sector: float(value) for sector, value in zip(sectors, lam)},
        "shares": {sector: float(value) for sector, value in zip(sectors, shares)},
        "targets": {sector: float(value) for sector, value in zip(sectors, targets)},
        "loss": fit_loss,
        "success": success,
        "best_fraction": best_fraction,
    }


def write_grid(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_growth(data: dict, cfg: Config, out: Path) -> None:
    z_plot = np.linspace(0.0, 0.8, 400)
    plt.figure(figsize=(9, 5.5))
    plt.plot(z_plot, fsigma8_ref(z_plot, cfg), lw=2.5, label="Planck LCDM proxy")
    plt.errorbar(
        data["z"],
        data["fsigma8_obs"],
        yerr=data["fsigma8_err"],
        fmt="o",
        capsize=3,
        label=r"BOSS DR12 $f\sigma_8$",
    )
    plt.xlabel("Redshift z")
    plt.ylabel(r"$f\sigma_8(z)$")
    plt.title(r"Growth-connectivity input: measured $f\sigma_8$")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=240)
    plt.close()


def plot_supports(z: np.ndarray, d: dict[str, np.ndarray], out: Path) -> None:
    plt.figure(figsize=(10, 6))
    plt.plot(z, d["Gamma_H"], label=r"$\Gamma_H$ expansion consistency")
    plt.plot(z, d["Gamma_V"], label=r"$V_{\rm corr}$ growth connectivity")
    plt.plot(z, d["Gamma_R"], label=r"$R_{\rm robust}$ joint robustness")
    plt.xlabel("Redshift z")
    plt.ylabel("Support")
    plt.ylim(0.0, 1.05)
    plt.title("Operational H, V, and R support sectors")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=240)
    plt.close()


def plot_cci(z: np.ndarray, d: dict[str, np.ndarray], zc: float, out: Path) -> None:
    plt.figure(figsize=(10, 6))
    plt.semilogy(z, d["CCI_v02_norm"], lw=2.5, label="v0.2 H-only projection")
    plt.semilogy(z, d["CCI_v03_norm"], lw=2.5, label="v0.3 H + V/R projection")
    plt.axvline(zc, color="black", ls="--", lw=1.5, label=rf"$z_c={zc:.3f}$")
    plt.xlabel("Redshift z")
    plt.ylabel("Normalised CCI")
    plt.title("Cosmological CCI: v0.2 vs v0.3")
    plt.grid(alpha=0.25, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=240)
    plt.close()


def plot_curvature(z: np.ndarray, curvature: np.ndarray, zc: float, out: Path) -> None:
    plt.figure(figsize=(10, 5))
    plt.plot(z, curvature, lw=2.0)
    plt.axhline(0.0, color="black", lw=1.0)
    plt.axvline(zc, color="red", ls="--", label=rf"$z_c={zc:.3f}$")
    plt.xlabel("Redshift z")
    plt.ylabel(r"$d^2\log \widehat{\mathrm{CCI}}_{\rm cos}/dz^2$")
    plt.title("Data-driven transition proxy from CCI curvature")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=240)
    plt.close()


def plot_lambdas(fit: dict, out: Path) -> None:
    sectors = fit["sectors"]
    values = [fit["lambda"][s] for s in sectors]
    plt.figure(figsize=(7, 4.8))
    plt.bar(sectors, values)
    plt.ylabel(r"$\lambda_a$")
    plt.title("MaxEnt-calibrated companion sector weights")
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(out, dpi=240)
    plt.close()


def main() -> None:
    cfg = Config()
    base = Path(__file__).resolve().parent
    plot_dir = base / "plots"
    plot_dir.mkdir(exist_ok=True)

    chron = read_csv_columns(base / "maat_cci_cosmology_v03_chronometers.csv")
    growth = read_csv_columns(base / "maat_cci_cosmology_v03_fsigma8.csv")

    z = np.linspace(cfg.z_min, cfg.z_max, cfg.n_grid)
    e = e_lcdm(z, cfg)
    d_growth = growth_factor(z, cfg)
    f8_ref_grid = fsigma8_ref(z, cfg)

    h_ref_data = cfg.H0 * e_lcdm(chron["z"], cfg)
    h_pull = (chron["H_obs_km_s_Mpc"] - h_ref_data) / chron["H_err"]
    chi2_h = local_chi2(z, chron["z"], h_pull, cfg.kernel_h)

    f8_ref_data = fsigma8_ref(growth["z"], cfg)
    f8_pull = (growth["fsigma8_obs"] - f8_ref_data) / growth["fsigma8_err"]
    chi2_v = local_chi2(z, growth["z"], f8_pull, cfg.kernel_growth)

    gamma_h = 1.0 / (1.0 + chi2_h)
    gamma_v = 1.0 / (1.0 + chi2_v)
    chi2_r = 0.5 * (chi2_h + chi2_v)
    gamma_r = 1.0 / (1.0 + chi2_r)

    u_struct = np.abs(e - d_growth) / (e + d_growth + EPS)
    cci2 = cci_v02(z, cfg)
    cci3 = cci_v03(z, cfg, gamma_v, gamma_r)
    cci2_norm = cci2 / cci2[0]
    cci3_norm = cci3 / cci3[0]

    zc, curvature = transition_redshift(z, cci3_norm)

    d_h = chi2_h
    d_s = u_struct
    d_v = chi2_v
    d_r = chi2_r
    defects = np.column_stack([d_h, d_s, d_v, d_r])
    lambda_fit = fit_maxent_lambdas(
        defects,
        cfg.maxent_best_fraction,
        cfg.maxent_prior_strength,
        cfg.maxent_balance_strength,
    )
    lam_vec = np.array([lambda_fit["lambda"][s] for s in lambda_fit["sectors"]])
    f_struct = defects @ lam_vec

    grid_rows = []
    for idx in range(len(z)):
        grid_rows.append(
            {
                "z": f"{z[idx]:.10g}",
                "E_LCDM": f"{e[idx]:.10g}",
                "D_growth": f"{d_growth[idx]:.10g}",
                "fsigma8_ref": f"{f8_ref_grid[idx]:.10g}",
                "U_struct": f"{u_struct[idx]:.10g}",
                "d_H": f"{d_h[idx]:.10g}",
                "d_S": f"{d_s[idx]:.10g}",
                "d_V": f"{d_v[idx]:.10g}",
                "d_R": f"{d_r[idx]:.10g}",
                "Gamma_H": f"{gamma_h[idx]:.10g}",
                "Gamma_V": f"{gamma_v[idx]:.10g}",
                "Gamma_R": f"{gamma_r[idx]:.10g}",
                "CCI_v02_norm": f"{cci2_norm[idx]:.10g}",
                "CCI_v03_norm": f"{cci3_norm[idx]:.10g}",
                "F_struct_v03": f"{f_struct[idx]:.10g}",
                "curvature_log_CCI_v03": f"{curvature[idx]:.10g}",
            }
        )
    write_grid(base / "maat_cci_cosmology_v03_grid.csv", grid_rows)

    payload = {
        "status": "Diagnostic v0.3 benchmark; not a precision cosmology likelihood.",
        "config": asdict(cfg),
        "transition_redshift_zc": zc,
        "lambda_fit": lambda_fit,
        "h_pull_rms": float(np.sqrt(np.mean(h_pull**2))),
        "fsigma8_pull_rms": float(np.sqrt(np.mean(f8_pull**2))),
        "cci_v03_norm_at_z1": float(np.interp(1.0, z, cci3_norm)),
        "cci_v03_norm_at_z2": float(np.interp(2.0, z, cci3_norm)),
        "min_R_support": float(np.min(gamma_r)),
        "min_V_support": float(np.min(gamma_v)),
    }
    (base / "maat_cci_cosmology_v03_results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    diagnostics = {
        "Gamma_H": gamma_h,
        "Gamma_V": gamma_v,
        "Gamma_R": gamma_r,
        "CCI_v02_norm": cci2_norm,
        "CCI_v03_norm": cci3_norm,
    }
    plot_growth(growth, cfg, plot_dir / "fsigma8_growth_connectivity.png")
    plot_supports(z, diagnostics, plot_dir / "v_r_supports.png")
    plot_cci(z, diagnostics, zc, plot_dir / "cci_v02_v03_comparison.png")
    plot_curvature(z, curvature, zc, plot_dir / "transition_curvature.png")
    plot_lambdas(lambda_fit, plot_dir / "lambda_calibration.png")

    print("Cosmological CCI v0.3 complete.")
    print(f"z_c = {zc:.4f}")
    print("lambda =", lambda_fit["lambda"])
    print("Wrote maat_cci_cosmology_v03_grid.csv")
    print("Wrote maat_cci_cosmology_v03_results.json")
    print("Wrote plots/*.png")


if __name__ == "__main__":
    main()

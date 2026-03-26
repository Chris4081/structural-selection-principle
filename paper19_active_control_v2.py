import os
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple


# ============================================================
# Paper 19 - Active Structural Control in 2D phi^4
# Version 2
# ------------------------------------------------------------
# Features:
# - 2D damped phi^4 field evolution
# - No-control vs control comparison
# - CCI proxy, structural free energy, mutual information
# - CSV export of time series
# - automatic figure saving
# - summary TXT + config JSON
# - final field snapshots as .npy
# ============================================================


# -----------------------------
# Configuration
# -----------------------------
@dataclass
class SimConfig:
    nx: int = 64
    ny: int = 64
    dx: float = 1.0
    dt: float = 0.02
    steps: int = 2500
    damping: float = 0.06
    seed: int = 42

    diag_every: int = 20
    hist_bins_mi: int = 24
    hist_bins_entropy: int = 32

    # CCI proxy parameters
    lambda_pi: float = 1.0
    s0: float = 1.0
    alpha_grad: float = 0.5
    i0: float = 0.5
    lambda_bound: float = 1.0
    phi_max: float = 3.0
    kappa_u: float = 1.0
    eps: float = 1e-8

    # Structural free energy proxy
    lam_stab: float = 1.0
    lam_conn: float = 1.0
    lam_dyn: float = 1.0
    s_star: float = 0.42
    sigma_s: float = 0.18

    # Control
    control_mode: str = "none"        # "none", "domainwall", "critical_cci"
    control_strength: float = 0.0
    target_cci: float = 0.32
    wall_width: float = 6.0

    # Naming
    run_name: str = "run"


# -----------------------------
# Utilities
# -----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def laplacian_periodic(phi: np.ndarray, dx: float) -> np.ndarray:
    return (
        np.roll(phi, +1, axis=0)
        + np.roll(phi, -1, axis=0)
        + np.roll(phi, +1, axis=1)
        + np.roll(phi, -1, axis=1)
        - 4.0 * phi
    ) / (dx * dx)


def grad_sq(phi: np.ndarray, dx: float) -> np.ndarray:
    dpx = (np.roll(phi, -1, axis=0) - np.roll(phi, +1, axis=0)) / (2.0 * dx)
    dpy = (np.roll(phi, -1, axis=1) - np.roll(phi, +1, axis=1)) / (2.0 * dx)
    return dpx**2 + dpy**2


def make_domainwall_target(nx: int, ny: int, width: float = 6.0) -> np.ndarray:
    x = np.arange(nx) - nx / 2.0
    X = np.repeat(x[:, None], ny, axis=1)
    return np.tanh(X / width).astype(np.float64)


def make_initial_field(
    nx: int,
    ny: int,
    rng: np.random.Generator,
    mode: str = "mixed",
) -> Tuple[np.ndarray, np.ndarray]:
    x = np.arange(nx) - nx / 2.0
    y = np.arange(ny) - ny / 2.0
    X, Y = np.meshgrid(x, y, indexing="ij")

    if mode == "chaotic":
        phi = rng.normal(0.0, 0.9, size=(nx, ny))
        pi = rng.normal(0.0, 0.4, size=(nx, ny))
        return phi, pi

    if mode == "domainwall":
        phi = np.tanh(X / 5.0)
        phi += 0.08 * rng.normal(size=(nx, ny))
        pi = 0.05 * rng.normal(size=(nx, ny))
        return phi, pi

    if mode == "localized":
        r2 = X**2 + Y**2
        phi = 1.4 * np.exp(-r2 / (2.0 * 7.0**2))
        phi += 0.04 * rng.normal(size=(nx, ny))
        pi = 0.03 * rng.normal(size=(nx, ny))
        return phi, pi

    # mixed
    a_vac = rng.uniform(0.0, 1.0)
    a_wall = rng.uniform(0.0, 1.0)
    a_loc = rng.uniform(0.0, 1.0)
    a_noise = rng.uniform(0.0, 1.0)

    vacuum = np.ones((nx, ny)) * rng.choice([-1.0, 1.0])
    wall = np.tanh(X / 5.0)
    localized = 1.2 * np.exp(-(X**2 + Y**2) / (2.0 * 7.0**2))
    noise = rng.normal(0.0, 1.0, size=(nx, ny))

    phi = a_vac * vacuum + a_wall * wall + a_loc * localized + 0.5 * a_noise * noise
    phi /= (1.0 + np.std(phi))
    phi = np.clip(phi, -2.5, 2.5)

    pi = 0.15 * rng.normal(size=(nx, ny))
    return phi, pi


# -----------------------------
# Entropy / information
# -----------------------------
def coarse_grained_entropy(phi: np.ndarray, bins: int = 32, phi_range=(-3, 3)) -> float:
    hist, _ = np.histogram(phi.ravel(), bins=bins, range=phi_range, density=False)
    p = hist.astype(np.float64)
    p /= max(p.sum(), 1.0)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p + 1e-12)))


def nearest_neighbor_mutual_information(
    phi: np.ndarray,
    bins: int = 24,
    phi_range=(-3, 3),
) -> float:
    a1 = phi
    b1 = np.roll(phi, -1, axis=0)

    a2 = phi
    b2 = np.roll(phi, -1, axis=1)

    a = np.concatenate([a1.ravel(), a2.ravel()])
    b = np.concatenate([b1.ravel(), b2.ravel()])

    joint, _, _ = np.histogram2d(
        a, b, bins=bins, range=[phi_range, phi_range], density=False
    )
    joint = joint.astype(np.float64)
    joint /= max(joint.sum(), 1.0)

    px = joint.sum(axis=1, keepdims=True)
    py = joint.sum(axis=0, keepdims=True)

    nz = joint > 0
    mi = np.sum(joint[nz] * np.log(joint[nz] / (px @ py)[nz] + 1e-12))
    return float(max(mi, 0.0))


# -----------------------------
# Diagnostics
# -----------------------------
def compute_diagnostics(
    phi: np.ndarray,
    pi: np.ndarray,
    cfg: SimConfig,
    s_prev: float = None,
) -> Dict[str, float]:
    lap = laplacian_periodic(phi, cfg.dx)
    g2 = grad_sq(phi, cfg.dx)

    pi2_mean = float(np.mean(pi**2))
    grad_mean = float(np.mean(g2))

    residual = lap - (phi**3 - phi)
    residual_abs_mean = float(np.mean(np.abs(residual)))

    s_now = coarse_grained_entropy(phi, bins=cfg.hist_bins_entropy, phi_range=(-3, 3))
    if s_prev is None:
        sdot_plus = 0.0
    else:
        sdot_plus = max(0.0, (s_now - s_prev) / cfg.dt)

    inn = nearest_neighbor_mutual_information(
        phi, bins=cfg.hist_bins_mi, phi_range=(-3, 3)
    )

    gamma_inst = pi2_mean / (cfg.lambda_pi + pi2_mean)
    s_dyn = (pi2_mean + cfg.alpha_grad * grad_mean) / (
        10.0 + pi2_mean + cfg.alpha_grad * grad_mean
    )
    gamma_prod = s_dyn
    gamma_coh = float(np.exp(-residual_abs_mean))
    gamma_corr = inn / (cfg.i0 + inn + cfg.eps)

    boundary_penalty = float(np.mean(np.maximum(0.0, np.abs(phi) - cfg.phi_max) ** 2))
    gamma_int = 1.0 - boundary_penalty / (cfg.lambda_bound + boundary_penalty + cfg.eps)
    gamma_int = float(np.clip(gamma_int, 0.0, 1.0))

    sectors = np.array([gamma_inst, gamma_prod, gamma_coh, gamma_corr, gamma_int])
    u_struct = float(np.std(sectors))

    cci = (
        gamma_inst
        * gamma_prod
        * (1.0 + cfg.kappa_u * u_struct)
        / (gamma_coh + gamma_corr + gamma_int + cfg.eps)
    )

    e_stab = gamma_inst + boundary_penalty / (cfg.lambda_bound + boundary_penalty + cfg.eps)
    e_conn = -inn
    e_dyn = ((s_dyn - cfg.s_star) / cfg.sigma_s) ** 2
    f_struct = cfg.lam_stab * e_stab + cfg.lam_conn * e_conn + cfg.lam_dyn * e_dyn

    return {
        "entropy": s_now,
        "sdot_plus": sdot_plus,
        "inn": inn,
        "pi2_mean": pi2_mean,
        "grad_mean": grad_mean,
        "residual_abs_mean": residual_abs_mean,
        "gamma_inst": gamma_inst,
        "gamma_prod": gamma_prod,
        "gamma_coh": gamma_coh,
        "gamma_corr": gamma_corr,
        "gamma_int": gamma_int,
        "u_struct": u_struct,
        "s_dyn": s_dyn,
        "cci": float(cci),
        "f_struct": float(f_struct),
        "boundary_penalty": boundary_penalty,
    }


# -----------------------------
# Control
# -----------------------------
def control_term(
    phi: np.ndarray,
    pi: np.ndarray,
    diag: Dict[str, float],
    cfg: SimConfig,
    target_field: np.ndarray = None,
) -> np.ndarray:
    if cfg.control_mode == "none" or cfg.control_strength == 0.0:
        return np.zeros_like(phi)

    if cfg.control_mode == "domainwall":
        if target_field is None:
            raise ValueError("target_field required for domainwall control mode.")
        return cfg.control_strength * (target_field - phi)

    if cfg.control_mode == "critical_cci":
        delta = cfg.target_cci - diag["cci"]
        local_shape = np.tanh(phi) - 0.5 * phi
        return cfg.control_strength * delta * local_shape

    raise ValueError(f"Unknown control_mode: {cfg.control_mode}")


# -----------------------------
# Dynamics
# -----------------------------
def rhs(
    phi: np.ndarray,
    pi: np.ndarray,
    cfg: SimConfig,
    diag: Dict[str, float],
    target_field: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray]:
    dphi = pi
    lap = laplacian_periodic(phi, cfg.dx)
    u = control_term(phi, pi, diag, cfg, target_field=target_field)
    dpi = lap - (phi**3 - phi) - cfg.damping * pi + u
    return dphi, dpi


def rk4_step(
    phi: np.ndarray,
    pi: np.ndarray,
    cfg: SimConfig,
    diag: Dict[str, float],
    target_field: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray]:
    dt = cfg.dt

    k1_phi, k1_pi = rhs(phi, pi, cfg, diag, target_field)

    phi2 = phi + 0.5 * dt * k1_phi
    pi2 = pi + 0.5 * dt * k1_pi
    diag2 = compute_diagnostics(phi2, pi2, cfg, s_prev=None)
    k2_phi, k2_pi = rhs(phi2, pi2, cfg, diag2, target_field)

    phi3 = phi + 0.5 * dt * k2_phi
    pi3 = pi + 0.5 * dt * k2_pi
    diag3 = compute_diagnostics(phi3, pi3, cfg, s_prev=None)
    k3_phi, k3_pi = rhs(phi3, pi3, cfg, diag3, target_field)

    phi4 = phi + dt * k3_phi
    pi4 = pi + dt * k3_pi
    diag4 = compute_diagnostics(phi4, pi4, cfg, s_prev=None)
    k4_phi, k4_pi = rhs(phi4, pi4, cfg, diag4, target_field)

    phi_new = phi + (dt / 6.0) * (k1_phi + 2 * k2_phi + 2 * k3_phi + k4_phi)
    pi_new = pi + (dt / 6.0) * (k1_pi + 2 * k2_pi + 2 * k3_pi + k4_pi)

    phi_new = np.clip(phi_new, -4.0, 4.0)
    pi_new = np.clip(pi_new, -4.0, 4.0)
    return phi_new, pi_new


# -----------------------------
# Run simulation
# -----------------------------
def run_simulation(cfg: SimConfig, init_mode: str = "chaotic") -> Dict[str, object]:
    rng = np.random.default_rng(cfg.seed)
    phi, pi = make_initial_field(cfg.nx, cfg.ny, rng, mode=init_mode)

    target_field = None
    if cfg.control_mode == "domainwall":
        target_field = make_domainwall_target(cfg.nx, cfg.ny, width=cfg.wall_width)

    t_hist: List[float] = []
    cci_hist: List[float] = []
    f_hist: List[float] = []
    inn_hist: List[float] = []
    sdot_hist: List[float] = []
    entropy_hist: List[float] = []
    pi2_hist: List[float] = []
    grad_hist: List[float] = []
    residual_hist: List[float] = []

    s_prev = None
    last_diag = None

    for step in range(cfg.steps):
        diag = compute_diagnostics(phi, pi, cfg, s_prev=s_prev)
        s_prev = diag["entropy"]

        phi, pi = rk4_step(phi, pi, cfg, diag, target_field=target_field)

        if step % cfg.diag_every == 0:
            last_diag = compute_diagnostics(phi, pi, cfg, s_prev=s_prev)
            t_hist.append(step * cfg.dt)
            cci_hist.append(last_diag["cci"])
            f_hist.append(last_diag["f_struct"])
            inn_hist.append(last_diag["inn"])
            sdot_hist.append(last_diag["sdot_plus"])
            entropy_hist.append(last_diag["entropy"])
            pi2_hist.append(last_diag["pi2_mean"])
            grad_hist.append(last_diag["grad_mean"])
            residual_hist.append(last_diag["residual_abs_mean"])

    if last_diag is None:
        last_diag = compute_diagnostics(phi, pi, cfg, s_prev=s_prev)

    return {
        "cfg": asdict(cfg),
        "phi_final": phi,
        "pi_final": pi,
        "target_field": target_field,
        "t": np.array(t_hist),
        "cci": np.array(cci_hist),
        "f_struct": np.array(f_hist),
        "inn": np.array(inn_hist),
        "sdot_plus": np.array(sdot_hist),
        "entropy": np.array(entropy_hist),
        "pi2_mean": np.array(pi2_hist),
        "grad_mean": np.array(grad_hist),
        "residual_abs_mean": np.array(residual_hist),
        "diag_final": last_diag,
    }


# -----------------------------
# Export helpers
# -----------------------------
def save_timeseries_csv(run: Dict[str, object], out_csv: str) -> None:
    rows = zip(
        run["t"],
        run["cci"],
        run["f_struct"],
        run["inn"],
        run["sdot_plus"],
        run["entropy"],
        run["pi2_mean"],
        run["grad_mean"],
        run["residual_abs_mean"],
    )
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "t",
            "cci",
            "f_struct",
            "inn",
            "sdot_plus",
            "entropy",
            "pi2_mean",
            "grad_mean",
            "residual_abs_mean",
        ])
        writer.writerows(rows)


def save_summary_txt(run: Dict[str, object], out_txt: str, label: str) -> None:
    d = run["diag_final"]
    lines = [
        f"=== {label} ===",
        f"Final CCI         : {d['cci']:.6f}",
        f"Final F_struct    : {d['f_struct']:.6f}",
        f"Final I_nn        : {d['inn']:.6f}",
        f"Final sdot_plus   : {d['sdot_plus']:.6f}",
        f"Final entropy     : {d['entropy']:.6f}",
        f"Final gamma_coh   : {d['gamma_coh']:.6f}",
        f"Final gamma_corr  : {d['gamma_corr']:.6f}",
        f"Final gamma_inst  : {d['gamma_inst']:.6f}",
        f"Final gamma_int   : {d['gamma_int']:.6f}",
        f"Final s_dyn       : {d['s_dyn']:.6f}",
        f"Residual mean     : {d['residual_abs_mean']:.6f}",
        f"Mean pi^2         : {d['pi2_mean']:.6f}",
        f"Mean grad^2       : {d['grad_mean']:.6f}",
        "",
        "Config:",
        json.dumps(run["cfg"], indent=2),
    ]
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def save_config_json(cfg: SimConfig, out_json: str) -> None:
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(asdict(cfg), f, indent=2)


def save_field_arrays(run: Dict[str, object], out_dir: str, prefix: str) -> None:
    np.save(os.path.join(out_dir, f"{prefix}_phi_final.npy"), run["phi_final"])
    np.save(os.path.join(out_dir, f"{prefix}_pi_final.npy"), run["pi_final"])
    if run.get("target_field", None) is not None:
        np.save(os.path.join(out_dir, f"{prefix}_target_field.npy"), run["target_field"])


# -----------------------------
# Plotting
# -----------------------------
def save_final_field_figure(
    run_no: Dict[str, object],
    run_ctrl: Dict[str, object],
    out_png: str,
    title_suffix: str = "",
) -> None:
    phi_a = run_no["phi_final"]
    phi_b = run_ctrl["phi_final"]
    target_b = run_ctrl.get("target_field", None)

    fig = plt.figure(figsize=(14, 4.8))

    ax1 = plt.subplot(1, 3, 1)
    im1 = ax1.imshow(phi_a, origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    ax1.set_title("Final field: no control")
    plt.colorbar(im1, ax=ax1, fraction=0.046)

    ax2 = plt.subplot(1, 3, 2)
    im2 = ax2.imshow(phi_b, origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    ax2.set_title("Final field: with control")
    plt.colorbar(im2, ax=ax2, fraction=0.046)

    ax3 = plt.subplot(1, 3, 3)
    if target_b is not None:
        im3 = ax3.imshow(target_b, origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
        ax3.set_title("Control target")
        plt.colorbar(im3, ax=ax3, fraction=0.046)
    else:
        ax3.axis("off")
        ax3.text(0.5, 0.5, "No target field", ha="center", va="center", fontsize=12)

    plt.suptitle(f"Paper 19 - Final Fields {title_suffix}", fontsize=14)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def save_timeseries_figure(
    run_no: Dict[str, object],
    run_ctrl: Dict[str, object],
    out_png: str,
    title_suffix: str = "",
) -> None:
    fig = plt.figure(figsize=(14, 10))

    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(run_no["t"], run_no["cci"], label="no control")
    ax1.plot(run_ctrl["t"], run_ctrl["cci"], label="control")
    ax1.set_title("CCI(t)")
    ax1.set_xlabel("t")
    ax1.set_ylabel("CCI")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(run_no["t"], run_no["f_struct"], label="no control")
    ax2.plot(run_ctrl["t"], run_ctrl["f_struct"], label="control")
    ax2.set_title("F_struct(t)")
    ax2.set_xlabel("t")
    ax2.set_ylabel("F_struct")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(run_no["t"], run_no["inn"], label="no control")
    ax3.plot(run_ctrl["t"], run_ctrl["inn"], label="control")
    ax3.set_title("Nearest-neighbor MI(t)")
    ax3.set_xlabel("t")
    ax3.set_ylabel("I_nn")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(run_no["t"], run_no["sdot_plus"], label="no control")
    ax4.plot(run_ctrl["t"], run_ctrl["sdot_plus"], label="control")
    ax4.set_title("Positive entropy-turnover proxy")
    ax4.set_xlabel("t")
    ax4.set_ylabel("sdot_plus")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.suptitle(f"Paper 19 - Structural Diagnostics {title_suffix}", fontsize=14)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def save_overlay_figure(
    run_no: Dict[str, object],
    run_ctrl: Dict[str, object],
    out_png: str,
    title_suffix: str = "",
) -> None:
    fig = plt.figure(figsize=(12, 8))

    plt.scatter(run_no["inn"], run_no["cci"], s=18, alpha=0.7, label="no control")
    plt.scatter(run_ctrl["inn"], run_ctrl["cci"], s=18, alpha=0.7, label="control")
    plt.xlabel("I_nn")
    plt.ylabel("CCI")
    plt.title(f"CCI vs I_nn {title_suffix}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def save_combined_comparison_csv(
    run_no: Dict[str, object],
    run_ctrl: Dict[str, object],
    out_csv: str,
) -> None:
    n = min(len(run_no["t"]), len(run_ctrl["t"]))
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "t",
            "cci_no",
            "cci_ctrl",
            "fstruct_no",
            "fstruct_ctrl",
            "inn_no",
            "inn_ctrl",
            "sdot_no",
            "sdot_ctrl",
        ])
        for i in range(n):
            writer.writerow([
                run_no["t"][i],
                run_no["cci"][i],
                run_ctrl["cci"][i],
                run_no["f_struct"][i],
                run_ctrl["f_struct"][i],
                run_no["inn"][i],
                run_ctrl["inn"][i],
                run_no["sdot_plus"][i],
                run_ctrl["sdot_plus"][i],
            ])


# -----------------------------
# Console summary
# -----------------------------
def print_summary(label: str, run: Dict[str, object]) -> None:
    d = run["diag_final"]
    print(f"\n=== {label} ===")
    print(f"Final CCI         : {d['cci']:.4f}")
    print(f"Final F_struct    : {d['f_struct']:.4f}")
    print(f"Final I_nn        : {d['inn']:.4f}")
    print(f"Final sdot_plus   : {d['sdot_plus']:.4f}")
    print(f"Final entropy     : {d['entropy']:.4f}")
    print(f"Final gamma_coh   : {d['gamma_coh']:.4f}")
    print(f"Final gamma_corr  : {d['gamma_corr']:.4f}")
    print(f"Final gamma_inst  : {d['gamma_inst']:.4f}")
    print(f"Final gamma_int   : {d['gamma_int']:.4f}")
    print(f"Final s_dyn       : {d['s_dyn']:.4f}")
    print(f"Residual mean     : {d['residual_abs_mean']:.4f}")
    print(f"Mean pi^2         : {d['pi2_mean']:.4f}")
    print(f"Mean grad^2       : {d['grad_mean']:.4f}")


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    project_dir = os.path.join("paper19_output", f"paper19_{timestamp}")
    figs_dir = os.path.join(project_dir, "figures")
    data_dir = os.path.join(project_dir, "data")
    summary_dir = os.path.join(project_dir, "summary")

    ensure_dir(project_dir)
    ensure_dir(figs_dir)
    ensure_dir(data_dir)
    ensure_dir(summary_dir)

    # Base config
    base = SimConfig(
        nx=64,
        ny=64,
        dx=1.0,
        dt=0.02,
        steps=2500,
        damping=0.06,
        seed=42,
        diag_every=20,
    )

    # No control
    cfg_no = SimConfig(**asdict(base))
    cfg_no.run_name = "no_control"
    cfg_no.control_mode = "none"
    cfg_no.control_strength = 0.0

    # Controlled
    cfg_ctrl = SimConfig(**asdict(base))
    cfg_ctrl.run_name = "domainwall_control"
    cfg_ctrl.control_mode = "domainwall"
    cfg_ctrl.control_strength = 0.18
    cfg_ctrl.wall_width = 6.0

    # Alternative example:
    # cfg_ctrl.run_name = "critical_cci_control"
    # cfg_ctrl.control_mode = "critical_cci"
    # cfg_ctrl.control_strength = 2.0
    # cfg_ctrl.target_cci = 0.32

    # Save configs
    save_config_json(cfg_no, os.path.join(summary_dir, "config_no_control.json"))
    save_config_json(cfg_ctrl, os.path.join(summary_dir, "config_control.json"))

    # Run simulations
    run_no = run_simulation(cfg_no, init_mode="chaotic")
    run_ctrl = run_simulation(cfg_ctrl, init_mode="chaotic")

    # Console output
    print_summary("NO CONTROL", run_no)
    print_summary("WITH CONTROL", run_ctrl)

    # Save per-run CSV
    save_timeseries_csv(run_no, os.path.join(data_dir, "timeseries_no_control.csv"))
    save_timeseries_csv(run_ctrl, os.path.join(data_dir, "timeseries_control.csv"))

    # Save combined CSV
    save_combined_comparison_csv(
        run_no,
        run_ctrl,
        os.path.join(data_dir, "timeseries_comparison.csv"),
    )

    # Save summaries
    save_summary_txt(run_no, os.path.join(summary_dir, "summary_no_control.txt"), "NO CONTROL")
    save_summary_txt(run_ctrl, os.path.join(summary_dir, "summary_control.txt"), "WITH CONTROL")

    # Save arrays
    save_field_arrays(run_no, data_dir, "no_control")
    save_field_arrays(run_ctrl, data_dir, "control")

    # Save figures
    title_suffix = f"({cfg_ctrl.control_mode})"
    save_final_field_figure(
        run_no,
        run_ctrl,
        os.path.join(figs_dir, "fig1_final_fields.png"),
        title_suffix=title_suffix,
    )
    save_timeseries_figure(
        run_no,
        run_ctrl,
        os.path.join(figs_dir, "fig2_diagnostics_timeseries.png"),
        title_suffix=title_suffix,
    )
    save_overlay_figure(
        run_no,
        run_ctrl,
        os.path.join(figs_dir, "fig3_cci_vs_inn.png"),
        title_suffix=title_suffix,
    )

    # Master summary
    master_summary = {
        "project_dir": project_dir,
        "figures": [
            "fig1_final_fields.png",
            "fig2_diagnostics_timeseries.png",
            "fig3_cci_vs_inn.png",
        ],
        "data_files": [
            "timeseries_no_control.csv",
            "timeseries_control.csv",
            "timeseries_comparison.csv",
            "no_control_phi_final.npy",
            "no_control_pi_final.npy",
            "control_phi_final.npy",
            "control_pi_final.npy",
        ],
        "control_mode": cfg_ctrl.control_mode,
        "control_strength": cfg_ctrl.control_strength,
        "seed": cfg_ctrl.seed,
    }

    with open(os.path.join(project_dir, "paper19_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(master_summary, f, indent=2)

    print("\nSaved Paper 19 output to:")
    print(project_dir)
    print("\nFigures:")
    print(os.path.join(figs_dir, "fig1_final_fields.png"))
    print(os.path.join(figs_dir, "fig2_diagnostics_timeseries.png"))
    print(os.path.join(figs_dir, "fig3_cci_vs_inn.png"))
    print("\nCSVs:")
    print(os.path.join(data_dir, "timeseries_no_control.csv"))
    print(os.path.join(data_dir, "timeseries_control.csv"))
    print(os.path.join(data_dir, "timeseries_comparison.csv"))

    
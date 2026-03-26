import os
import csv
import json
import time
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict, replace
from typing import Dict, List, Tuple


# ============================================================
# Paper 19 - Active Structural Control in 2D phi^4
# Version 3: Parameter sweep over damping + control_strength
# ------------------------------------------------------------
# Features:
# - 2D damped phi^4 field evolution
# - sweep over damping and control strength
# - no-control baseline for each damping
# - controlled run for each (damping, control_strength)
# - exports:
#     * CSV results table
#     * heatmaps
#     * best-run summary
#     * best-run figures
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
    steps: int = 1800
    damping: float = 0.06
    seed: int = 42

    diag_every: int = 20
    hist_bins_mi: int = 24
    hist_bins_entropy: int = 32

    # CCI proxy parameters
    lambda_pi: float = 1.0
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
    mode: str = "chaotic",
) -> Tuple[np.ndarray, np.ndarray]:
    x = np.arange(nx) - nx / 2.0
    y = np.arange(ny) - ny / 2.0
    X, Y = np.meshgrid(x, y, indexing="ij")

    if mode == "chaotic":
        phi = rng.normal(0.0, 1.2, size=(nx, ny))
        pi = rng.normal(0.0, 0.7, size=(nx, ny))
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

    phi = rng.normal(0.0, 0.9, size=(nx, ny))
    pi = rng.normal(0.0, 0.3, size=(nx, ny))
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
# Simulation runner
# -----------------------------
def run_simulation(cfg: SimConfig, init_mode: str = "chaotic") -> Dict[str, object]:
    rng = np.random.default_rng(cfg.seed)
    phi, pi = make_initial_field(cfg.nx, cfg.ny, rng, mode=init_mode)

    target_field = None
    if cfg.control_mode == "domainwall":
        target_field = make_domainwall_target(cfg.nx, cfg.ny, width=cfg.wall_width)

    t_hist = []
    cci_hist = []
    f_hist = []
    inn_hist = []

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
        "diag_final": last_diag,
    }


# -----------------------------
# Plot helpers
# -----------------------------
def save_heatmap(
    matrix: np.ndarray,
    xvals: List[float],
    yvals: List[float],
    xlabel: str,
    ylabel: str,
    title: str,
    out_png: str,
) -> None:
    fig = plt.figure(figsize=(8, 6))
    plt.imshow(matrix, origin="lower", aspect="auto")
    plt.colorbar()
    plt.xticks(range(len(xvals)), [f"{x:.2f}" for x in xvals], rotation=45)
    plt.yticks(range(len(yvals)), [f"{y:.2f}" for y in yvals])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def save_best_run_figure(
    baseline_run: Dict[str, object],
    best_run: Dict[str, object],
    out_png: str,
    title_suffix: str = "",
) -> None:
    fig = plt.figure(figsize=(14, 10))

    ax1 = plt.subplot(2, 3, 1)
    im1 = ax1.imshow(baseline_run["phi_final"], origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    ax1.set_title("Baseline final field")
    plt.colorbar(im1, ax=ax1, fraction=0.046)

    ax2 = plt.subplot(2, 3, 2)
    im2 = ax2.imshow(best_run["phi_final"], origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
    ax2.set_title("Best controlled field")
    plt.colorbar(im2, ax=ax2, fraction=0.046)

    ax3 = plt.subplot(2, 3, 3)
    if best_run.get("target_field", None) is not None:
        im3 = ax3.imshow(best_run["target_field"], origin="lower", cmap="coolwarm", vmin=-1.5, vmax=1.5)
        ax3.set_title("Target field")
        plt.colorbar(im3, ax=ax3, fraction=0.046)
    else:
        ax3.axis("off")

    ax4 = plt.subplot(2, 3, 4)
    ax4.plot(baseline_run["t"], baseline_run["cci"], label="baseline")
    ax4.plot(best_run["t"], best_run["cci"], label="best controlled")
    ax4.set_title("CCI(t)")
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    ax5 = plt.subplot(2, 3, 5)
    ax5.plot(baseline_run["t"], baseline_run["f_struct"], label="baseline")
    ax5.plot(best_run["t"], best_run["f_struct"], label="best controlled")
    ax5.set_title("F_struct(t)")
    ax5.grid(True, alpha=0.3)
    ax5.legend()

    ax6 = plt.subplot(2, 3, 6)
    ax6.plot(baseline_run["t"], baseline_run["inn"], label="baseline")
    ax6.plot(best_run["t"], best_run["inn"], label="best controlled")
    ax6.set_title("I_nn(t)")
    ax6.grid(True, alpha=0.3)
    ax6.legend()

    plt.suptitle(f"Best Paper 19 Control Run {title_suffix}", fontsize=15)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Scoring
# -----------------------------
def compute_score(
    delta_cci: float,
    delta_f: float,
    delta_inn: float,
    baseline_cci: float,
    controlled_cci: float,
) -> float:
    # We want:
    # - lower CCI
    # - lower F_struct
    # - higher I_nn
    # Extra weight if baseline CCI is not already ultra-tiny
    regime_weight = 1.0 + min(1.0, 5.0 * baseline_cci)
    cci_term = -delta_cci
    f_term = -0.5 * delta_f
    inn_term = 1.5 * delta_inn
    return regime_weight * (cci_term + f_term + inn_term)


# -----------------------------
# Main sweep
# -----------------------------
if __name__ == "__main__":
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    project_dir = os.path.join("paper19_output", f"paper19_sweep_{timestamp}")
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
        steps=1800,
        seed=42,
        diag_every=20,
        control_mode="domainwall",
        wall_width=6.0,
    )

    # Sweep values
    damping_values = [0.02, 0.03, 0.04, 0.06, 0.08]
    control_values = [0.00, 0.08, 0.16, 0.24, 0.32]

    with open(os.path.join(summary_dir, "sweep_config.json"), "w", encoding="utf-8") as f:
        json.dump(
            {
                "base_config": asdict(base),
                "damping_values": damping_values,
                "control_values": control_values,
            },
            f,
            indent=2,
        )

    results = []
    baseline_runs_by_damping = {}

    print("\n=== Running baselines by damping ===")
    for damping in damping_values:
        cfg_base = replace(
            base,
            damping=damping,
            control_mode="none",
            control_strength=0.0,
            run_name=f"baseline_damp_{damping:.2f}",
        )
        run_base = run_simulation(cfg_base, init_mode="chaotic")
        baseline_runs_by_damping[damping] = run_base
        d = run_base["diag_final"]
        print(
            f"damping={damping:.2f} | baseline final CCI={d['cci']:.4f}, "
            f"F_struct={d['f_struct']:.4f}, I_nn={d['inn']:.4f}"
        )

    print("\n=== Running controlled sweep ===")
    best_score = -1e18
    best_entry = None
    best_run = None

    for iy, damping in enumerate(damping_values):
        baseline_run = baseline_runs_by_damping[damping]
        baseline_diag = baseline_run["diag_final"]

        for ix, control_strength in enumerate(control_values):
            if control_strength == 0.0:
                run_ctrl = baseline_run
            else:
                cfg_ctrl = replace(
                    base,
                    damping=damping,
                    control_mode="domainwall",
                    control_strength=control_strength,
                    run_name=f"ctrl_damp_{damping:.2f}_lam_{control_strength:.2f}",
                )
                run_ctrl = run_simulation(cfg_ctrl, init_mode="chaotic")

            ctrl_diag = run_ctrl["diag_final"]

            delta_cci = ctrl_diag["cci"] - baseline_diag["cci"]
            delta_f = ctrl_diag["f_struct"] - baseline_diag["f_struct"]
            delta_inn = ctrl_diag["inn"] - baseline_diag["inn"]

            score = compute_score(
                delta_cci=delta_cci,
                delta_f=delta_f,
                delta_inn=delta_inn,
                baseline_cci=baseline_diag["cci"],
                controlled_cci=ctrl_diag["cci"],
            )

            entry = {
                "damping": damping,
                "control_strength": control_strength,
                "baseline_cci": baseline_diag["cci"],
                "baseline_f_struct": baseline_diag["f_struct"],
                "baseline_inn": baseline_diag["inn"],
                "controlled_cci": ctrl_diag["cci"],
                "controlled_f_struct": ctrl_diag["f_struct"],
                "controlled_inn": ctrl_diag["inn"],
                "delta_cci": delta_cci,
                "delta_f_struct": delta_f,
                "delta_inn": delta_inn,
                "score": score,
            }
            results.append(entry)

            print(
                f"damp={damping:.2f}, lambda={control_strength:.2f} | "
                f"dCCI={delta_cci:+.4e}, dF={delta_f:+.4f}, dInn={delta_inn:+.4f}, score={score:+.4f}"
            )

            if control_strength > 0.0 and score > best_score:
                best_score = score
                best_entry = entry
                best_run = run_ctrl

    # Save results CSV
    csv_path = os.path.join(data_dir, "paper19_sweep_results.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    # Build heatmaps
    ny = len(damping_values)
    nx = len(control_values)

    heat_delta_cci = np.zeros((ny, nx))
    heat_delta_f = np.zeros((ny, nx))
    heat_delta_inn = np.zeros((ny, nx))
    heat_score = np.zeros((ny, nx))
    heat_final_cci = np.zeros((ny, nx))

    for entry in results:
        iy = damping_values.index(entry["damping"])
        ix = control_values.index(entry["control_strength"])
        heat_delta_cci[iy, ix] = entry["delta_cci"]
        heat_delta_f[iy, ix] = entry["delta_f_struct"]
        heat_delta_inn[iy, ix] = entry["delta_inn"]
        heat_score[iy, ix] = entry["score"]
        heat_final_cci[iy, ix] = entry["controlled_cci"]

    save_heatmap(
        heat_delta_cci,
        control_values,
        damping_values,
        xlabel="control_strength",
        ylabel="damping",
        title="Delta CCI (controlled - baseline)",
        out_png=os.path.join(figs_dir, "heatmap_delta_cci.png"),
    )
    save_heatmap(
        heat_delta_f,
        control_values,
        damping_values,
        xlabel="control_strength",
        ylabel="damping",
        title="Delta F_struct (controlled - baseline)",
        out_png=os.path.join(figs_dir, "heatmap_delta_fstruct.png"),
    )
    save_heatmap(
        heat_delta_inn,
        control_values,
        damping_values,
        xlabel="control_strength",
        ylabel="damping",
        title="Delta I_nn (controlled - baseline)",
        out_png=os.path.join(figs_dir, "heatmap_delta_inn.png"),
    )
    save_heatmap(
        heat_score,
        control_values,
        damping_values,
        xlabel="control_strength",
        ylabel="damping",
        title="Composite control score",
        out_png=os.path.join(figs_dir, "heatmap_score.png"),
    )
    save_heatmap(
        heat_final_cci,
        control_values,
        damping_values,
        xlabel="control_strength",
        ylabel="damping",
        title="Final controlled CCI",
        out_png=os.path.join(figs_dir, "heatmap_final_cci.png"),
    )

    # Save best run
    if best_entry is not None:
        baseline_best = baseline_runs_by_damping[best_entry["damping"]]

        np.save(os.path.join(data_dir, "best_phi_final.npy"), best_run["phi_final"])
        np.save(os.path.join(data_dir, "best_pi_final.npy"), best_run["pi_final"])
        if best_run.get("target_field", None) is not None:
            np.save(os.path.join(data_dir, "best_target_field.npy"), best_run["target_field"])

        save_best_run_figure(
            baseline_best,
            best_run,
            os.path.join(figs_dir, "best_run_comparison.png"),
            title_suffix=f"(damping={best_entry['damping']:.2f}, lambda={best_entry['control_strength']:.2f})",
        )

        with open(os.path.join(summary_dir, "best_run_summary.txt"), "w", encoding="utf-8") as f:
            f.write("=== BEST PAPER 19 SWEEP RESULT ===\n")
            for k, v in best_entry.items():
                if isinstance(v, float):
                    f.write(f"{k}: {v:.8f}\n")
                else:
                    f.write(f"{k}: {v}\n")

    manifest = {
        "project_dir": project_dir,
        "csv_results": "data/paper19_sweep_results.csv",
        "figures": [
            "figures/heatmap_delta_cci.png",
            "figures/heatmap_delta_fstruct.png",
            "figures/heatmap_delta_inn.png",
            "figures/heatmap_score.png",
            "figures/heatmap_final_cci.png",
            "figures/best_run_comparison.png",
        ],
        "best_entry": best_entry,
    }

    with open(os.path.join(project_dir, "paper19_sweep_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print("\n=== Sweep finished ===")
    print("Output folder:")
    print(project_dir)
    print("\nBest entry:")
    print(best_entry)
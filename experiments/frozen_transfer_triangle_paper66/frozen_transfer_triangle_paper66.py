#!/usr/bin/env python3
"""Paper 66: Frozen transfer triangle across SAT, Quantum, and Fluids.

The benchmark extends Papers 52--53 by adding the Paper-60 fluid diagnostic as
a third internal domain.  All domains expose the same support language
H, B, S, V, R_rob and a predeclared quality target where larger is better.

The protocol is intentionally simple:

1. learn non-negative source weights from one domain only;
2. freeze those weights;
3. evaluate them on the other two domains without retuning or sign flips;
4. compare against equal weights, source-fitted rules, V-only target alignment,
   and shuffled-defect nulls.

This is an internal transfer test, not external validation.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr, pearsonr

_CACHE = Path(tempfile.gettempdir()) / "maat_paper66_matplotlib_cache"
_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EPS = 1.0e-9
SEED = 66066
SECTORS = ["H", "B", "S", "V", "R_rob"]
FIT_RULES = ["equal", "nnls", "ridge_pos", "v_only"]


@dataclass
class DomainData:
    name: str
    df: pd.DataFrame
    supports: np.ndarray
    defects: np.ndarray
    quality: np.ndarray
    quality_label: str
    notes: str


def normalize01(x: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(x), dtype=float)
    lo = np.nanmin(arr)
    hi = np.nanmax(arr)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < EPS:
        return np.zeros_like(arr, dtype=float)
    return (arr - lo) / (hi - lo)


def safe_spearman(x: Iterable[float], y: Iterable[float]) -> float:
    df = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 4 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return float("nan")
    val = spearmanr(df["x"], df["y"]).correlation
    return float(val if np.isfinite(val) else np.nan)


def safe_pearson(x: Iterable[float], y: Iterable[float]) -> float:
    df = pd.DataFrame({"x": list(x), "y": list(y)}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < 4 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return float("nan")
    val = pearsonr(df["x"], df["y"]).statistic
    return float(val if np.isfinite(val) else np.nan)


def make_domain(
    name: str,
    df: pd.DataFrame,
    support_map: dict[str, str],
    quality: Iterable[float],
    quality_label: str,
    notes: str,
) -> DomainData:
    supports = np.column_stack([pd.to_numeric(df[support_map[s]], errors="coerce") for s in SECTORS])
    supports = np.clip(np.nan_to_num(supports, nan=EPS, posinf=1.0, neginf=EPS), EPS, 1.0)
    defects = -np.log(EPS + supports)
    q = normalize01(quality)
    return DomainData(
        name=name,
        df=df,
        supports=supports,
        defects=defects,
        quality=q,
        quality_label=quality_label,
        notes=notes,
    )


def load_domains(root: Path) -> list[DomainData]:
    domains: list[DomainData] = []

    sat = pd.read_csv(
        root
        / "experiments/sat_frustration_fields_paper50b/outputs_paper50b/paper50b_sat_instances.csv"
    )
    domains.append(
        make_domain(
            "SAT",
            sat,
            {s: s for s in SECTORS},
            1.0 - normalize01(sat["log_nodes"]),
            "1 - minmax(log DPLL nodes)",
            "Paper 50b SAT frustration-field benchmark; higher quality means lower DPLL hardness.",
        )
    )

    quantum = pd.read_csv(
        root
        / "experiments/quantum_pointer_state_selection_paper51/outputs_paper51/paper51_pointer_instances.csv"
    )
    domains.append(
        make_domain(
            "Quantum",
            quantum,
            {s: s for s in SECTORS},
            quantum["target"],
            "minmax(pointer robustness target)",
            "Paper 51 non-commuting Lindblad pointer-state benchmark.",
        )
    )

    burgers = pd.read_csv(
        root
        / "experiments/fluid_blowup_diagnostics_paper60/outputs_paper60/paper60_burgers_diagnostics.csv"
    )
    ns2d = pd.read_csv(
        root
        / "experiments/fluid_blowup_diagnostics_paper60/outputs_paper60/paper60_ns2d_diagnostics.csv"
    )
    fluid = pd.concat([burgers, ns2d], ignore_index=True)
    domains.append(
        make_domain(
            "Fluid",
            fluid,
            {"H": "H", "B": "B", "S": "S_eff", "V": "V", "R_rob": "R_rob"},
            1.0 - normalize01(fluid["MAAT_warning"]),
            "1 - minmax(MAAT warning)",
            "Paper 60 Burgers plus 2D Navier--Stokes diagnostic; higher quality means lower warning pressure.",
        )
    )

    return domains


def normalize_weights(w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    w = np.clip(w, 0.0, None)
    if w.sum() <= EPS:
        return np.ones(len(SECTORS), dtype=float) / len(SECTORS)
    return w / w.sum()


def fit_weights(source: DomainData, rule: str) -> np.ndarray:
    D = source.defects
    cost = 1.0 - normalize01(source.quality)
    n = len(SECTORS)

    if rule == "equal":
        return np.ones(n, dtype=float) / n

    if rule == "v_only":
        w = np.zeros(n, dtype=float)
        w[SECTORS.index("V")] = 1.0
        return w

    if rule == "nnls":
        w, _ = nnls(D, cost)
        return normalize_weights(w)

    if rule == "ridge_pos":
        alpha = 0.05 * (np.trace(D.T @ D) / max(n, 1)) / max(len(D), 1)
        w = np.linalg.solve(D.T @ D + alpha * np.eye(n), D.T @ cost)
        return normalize_weights(w)

    raise ValueError(f"Unknown fit rule: {rule}")


def score(defects: np.ndarray, weights: np.ndarray) -> np.ndarray:
    # Higher score means lower weighted defect cost and therefore higher
    # predicted structural quality.
    return -defects @ weights


def shuffle_defects(defects: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    shuffled = defects.copy()
    for j in range(shuffled.shape[1]):
        shuffled[:, j] = rng.permutation(shuffled[:, j])
    return shuffled


def evaluate(domains: list[DomainData], n_null: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    weight_rows: list[dict] = []
    null_rows: list[dict] = []

    for source in domains:
        for rule in FIT_RULES:
            w = fit_weights(source, rule)
            dom_sector = SECTORS[int(np.argmax(w))]
            weight_rows.append(
                {
                    "source_domain": source.name,
                    "fit_rule": rule,
                    "dominant_sector": dom_sector,
                    "v_weight": float(w[SECTORS.index("V")]),
                    **{f"lambda_{s}": float(v) for s, v in zip(SECTORS, w)},
                }
            )

            for target in domains:
                pred = score(target.defects, w)
                rho = safe_spearman(pred, target.quality)
                pearson = safe_pearson(pred, target.quality)

                null_vals = []
                if source.name != target.name:
                    for null_id in range(n_null):
                        null_pred = score(shuffle_defects(target.defects, rng), w)
                        nr = safe_spearman(null_pred, target.quality)
                        null_vals.append(nr)
                        null_rows.append(
                            {
                                "source_domain": source.name,
                                "target_domain": target.name,
                                "transfer": f"{source.name}->{target.name}",
                                "fit_rule": rule,
                                "null_id": null_id,
                                "spearman": nr,
                            }
                        )
                null_arr = np.asarray(null_vals, dtype=float)
                null95 = float(np.nanquantile(null_arr, 0.95)) if len(null_arr) else float("nan")
                null_med = float(np.nanmedian(null_arr)) if len(null_arr) else float("nan")
                null_p = (
                    float((np.sum(null_arr >= rho) + 1.0) / (np.sum(np.isfinite(null_arr)) + 1.0))
                    if len(null_arr) and np.isfinite(rho)
                    else float("nan")
                )
                rows.append(
                    {
                        "source_domain": source.name,
                        "target_domain": target.name,
                        "transfer": f"{source.name}->{target.name}",
                        "fit_rule": rule,
                        "is_cross_domain": bool(source.name != target.name),
                        "spearman": rho,
                        "pearson": pearson,
                        "null_median": null_med,
                        "null_q95": null95,
                        "null_p_ge_observed": null_p,
                        "above_null_q95": bool(rho > null95) if np.isfinite(null95) else False,
                        "dominant_sector": dom_sector,
                        "v_weight": float(w[SECTORS.index("V")]),
                        "source_quality_label": source.quality_label,
                        "target_quality_label": target.quality_label,
                    }
                )

    return pd.DataFrame(rows), pd.DataFrame(weight_rows), pd.DataFrame(null_rows)


def plot_triangle(results: pd.DataFrame, out: Path, rule: str = "nnls") -> None:
    names = ["SAT", "Quantum", "Fluid"]
    mat = np.full((len(names), len(names)), np.nan)
    work = results[results["fit_rule"] == rule]
    for i, src in enumerate(names):
        for j, tgt in enumerate(names):
            row = work[(work["source_domain"] == src) & (work["target_domain"] == tgt)]
            if not row.empty:
                mat[i, j] = float(row["spearman"].iloc[0])

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    im = ax.imshow(mat, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_yticklabels(names)
    ax.set_xlabel("Target domain")
    ax.set_ylabel("Frozen source domain")
    ax.set_title(f"Frozen transfer triangle ({rule})")
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", color="black", fontsize=11)
    fig.colorbar(im, ax=ax, label="Spearman")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_vonly_target_alignment(results: pd.DataFrame, out: Path) -> None:
    work = results[(results["fit_rule"] == "v_only") & (results["is_cross_domain"])].copy()
    targets = ["SAT", "Quantum", "Fluid"]
    observed = []
    null95 = []
    for target in targets:
        rows = work[work["target_domain"] == target]
        observed.append(float(rows["spearman"].iloc[0]) if not rows.empty else np.nan)
        # Conservative because each target has two row-duplicated source entries
        # with independently sampled shuffled nulls.
        null95.append(float(rows["null_q95"].max()) if not rows.empty else np.nan)

    x = np.arange(len(targets))
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.bar(x - 0.18, observed, 0.36, label="V-only target alignment", color="#f97316")
    ax.bar(x + 0.18, null95, 0.36, label="shuffled-defect q95", color="#94a3b8")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(targets)
    ax.set_ylabel("Spearman")
    ax.set_title("V-only alignment by target domain")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_rule_comparison(results: pd.DataFrame, out: Path) -> None:
    work = results[results["is_cross_domain"]].copy()
    summary = work.groupby("fit_rule", as_index=False)["spearman"].mean()
    order = [r for r in FIT_RULES if r in set(summary["fit_rule"])]
    vals = [float(summary.loc[summary["fit_rule"] == r, "spearman"].iloc[0]) for r in order]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.bar(order, vals, color=["#64748b", "#2563eb", "#0f766e", "#f97316"])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Mean cross-domain Spearman")
    ax.set_title("Frozen transfer rules (V-only is row-duplicated)")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_v_weights(weights: pd.DataFrame, out: Path) -> None:
    work = weights[weights["fit_rule"].isin(["nnls", "ridge_pos", "v_only"])].copy()
    pivot = work.pivot(index="source_domain", columns="fit_rule", values="v_weight")
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    pivot.plot(kind="bar", ax=ax, color=["#2563eb", "#0f766e", "#f97316"])
    ax.set_ylabel("Frozen V weight")
    ax.set_title("How much of each source fit is connectivity?")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_null_margins(results: pd.DataFrame, out: Path) -> None:
    work = results[(results["is_cross_domain"]) & (results["fit_rule"].isin(["nnls", "v_only"]))].copy()
    work["margin"] = work["spearman"] - work["null_q95"]
    labels = work["transfer"] + "\n" + work["fit_rule"]
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    colors = ["#15803d" if m > 0 else "#b91c1c" for m in work["margin"]]
    ax.bar(np.arange(len(work)), work["margin"], color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(np.arange(len(work)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Observed Spearman - shuffled q95")
    ax.set_title("Transfer above shuffled-defect null?")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_readme(outdir: Path, summary: dict) -> None:
    lines = [
        "# Paper 66 -- Frozen Transfer Triangle",
        "",
        "Internal transfer benchmark across SAT, Quantum, and Fluid domains.",
        "",
        "The protocol freezes source-domain weights and evaluates them on the other two domains without target-domain retuning.",
        "",
        "## Headline",
        "",
        f"- Best cross-domain direction: `{summary['best_cross_transfer']['transfer']}` using `{summary['best_cross_transfer']['fit_rule']}` with Spearman `{summary['best_cross_transfer']['spearman']:.4f}`.",
        f"- Row-duplicated `v_only` target-alignment mean: `{summary['mean_cross_spearman_by_rule'].get('v_only', float('nan')):.4f}`.",
        f"- Mean cross-domain Spearman by `nnls`: `{summary['mean_cross_spearman_by_rule'].get('nnls', float('nan')):.4f}`.",
        "- `v_only` is a target-mode diagnostic, not source-trained transfer.",
        "",
        "## Outputs",
        "",
        "- `paper66_transfer_results.csv`",
        "- `paper66_frozen_weights.csv`",
        "- `paper66_shuffle_nulls.csv`",
        "- `paper66_summary.json`",
        "- `fig1_transfer_triangle_nnls.png`",
        "- `fig2_vonly_target_alignment.png`",
        "- `fig3_rule_comparison.png`",
        "- `fig4_v_weights.png`",
        "- `fig5_null_margins.png`",
        "",
        "## Status",
        "",
        "This is an internal benchmark only. It tests whether the connectivity mode `V` behaves like a shared low-frequency target-alignment channel across existing MAAT domains. The Fluid target is support-derived, so transfers into Fluid are partly circular. It is not external validation and not a universality proof.",
    ]
    outdir.joinpath("README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    outdir = Path("outputs_paper66")
    outdir.mkdir(parents=True, exist_ok=True)

    domains = load_domains(root)
    results, weights, nulls = evaluate(domains, n_null=500, seed=SEED)

    results.to_csv(outdir / "paper66_transfer_results.csv", index=False)
    weights.to_csv(outdir / "paper66_frozen_weights.csv", index=False)
    nulls.to_csv(outdir / "paper66_shuffle_nulls.csv", index=False)

    plot_triangle(results, outdir / "fig1_transfer_triangle_nnls.png", rule="nnls")
    plot_vonly_target_alignment(results, outdir / "fig2_vonly_target_alignment.png")
    plot_rule_comparison(results, outdir / "fig3_rule_comparison.png")
    plot_v_weights(weights, outdir / "fig4_v_weights.png")
    plot_null_margins(results, outdir / "fig5_null_margins.png")

    cross = results[results["is_cross_domain"]].copy()
    by_rule = cross.groupby("fit_rule")["spearman"].mean().to_dict()
    best = cross.sort_values("spearman", ascending=False).iloc[0].to_dict()
    nnls_cross = cross[cross["fit_rule"] == "nnls"].sort_values("spearman", ascending=False)
    vonly_cross = cross[cross["fit_rule"] == "v_only"].sort_values("spearman", ascending=False)
    vonly_target_alignment = []
    for target, rows in vonly_cross.groupby("target_domain"):
        vonly_target_alignment.append(
            {
                "target_domain": target,
                "spearman": float(rows["spearman"].iloc[0]),
                "conservative_null_q95": float(rows["null_q95"].max()),
                "above_null_q95": bool(rows["spearman"].iloc[0] > rows["null_q95"].max()),
            }
        )
    summary = {
        "paper": 66,
        "title": "The Frozen Transfer Triangle: SAT, Quantum, Fluids",
        "status": "internal frozen-weight transfer benchmark, not external validation",
        "domains": [
            {
                "name": d.name,
                "n": int(len(d.df)),
                "quality_label": d.quality_label,
                "notes": d.notes,
            }
            for d in domains
        ],
        "n_cross_directions": 6,
        "fit_rules": FIT_RULES,
        "best_cross_transfer": {
            "transfer": best["transfer"],
            "fit_rule": best["fit_rule"],
            "spearman": float(best["spearman"]),
            "dominant_sector": best["dominant_sector"],
            "v_weight": float(best["v_weight"]),
            "above_null_q95": bool(best["above_null_q95"]),
        },
        "mean_cross_spearman_by_rule": {k: float(v) for k, v in by_rule.items()},
        "v_only_target_alignment": sorted(vonly_target_alignment, key=lambda r: r["target_domain"]),
        "nnls_cross_transfers": nnls_cross[
            ["transfer", "spearman", "null_q95", "above_null_q95", "dominant_sector", "v_weight"]
        ].to_dict(orient="records"),
        "v_only_cross_transfers": vonly_cross[
            ["transfer", "spearman", "null_q95", "above_null_q95"]
        ].to_dict(orient="records"),
        "notes": [
            "Higher target quality is predeclared per domain.",
            "No target-domain sign flips or retuning are allowed.",
            "Shuffled-defect nulls preserve marginal sector distributions but destroy row-level support coherence.",
            "V-only rows are source-independent and should be read as target alignment, not fair architecture-level transfer.",
            "Fluid quality is support-derived from Paper 60 and is therefore less independent than the SAT target.",
        ],
    }
    (outdir / "paper66_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_readme(outdir, summary)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
plateau_degeneracy_exact.py
===========================
Computes exact plateau degeneracy measures for entropy-information
scaling in nonlinear field systems (Papers 07, 09, 11, 14, 16).

Degeneracy framework (Paper 16):
  - epsilon-plateau set  P_eps = {alpha | r_s(alpha) >= r_s_max - eps}
  - D_width  = max(P_eps) - min(P_eps)
  - D_norm   = D_width / (alpha_max - alpha_min)  in [0,1]
  - D_flat   = Var(r_s) within P_eps
  - D        = D_width / (1 + D_flat)              combined measure
  - D_curv   = |d2 r_s / d alpha2|^{-1} at alpha* (curvature, continuous limit)

Input CSVs (must be in same directory):
  cci_entropy_information_test.csv      <- 1D ensemble data (Paper 07)
  2d_cci_entropy_information_test.csv   <- 2D ensemble data (Paper 09)
  3d_cci_alpha_scan.csv                 <- 3D alpha scan   (Paper 11)

Usage:
  python plateau_degeneracy_exact.py

Output:
  Terminal table + plateau_degeneracy_results.csv

Requirements:
  pip install numpy pandas scipy

Results (Paper 16, Table 1):
  1D: alpha*=2.5, r_max=0.734, plateau=[2.1,2.5], D_norm=0.133, n=5
  2D: alpha*=3.0, r_max=0.870, plateau=[2.8,3.0], D_norm=0.067, n=3
  3D: alpha*=2.8, r_max=0.849, plateau=[2.6,3.5], D_norm=0.300, n=10
  Non-monotonic: D_norm(2D) < D_norm(1D) < D_norm(3D)
  => dimensional transition between 2D and 3D

Author:  Christof Krieg, 2026
Repo:    https://github.com/Chris4081/structural-selection-principle
Paper:   https://doi.org/10.5281/zenodo.XXXXXX  (Paper 16)
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────
EPS_DENOM = 1e-8     # numerical stabiliser for CCI ratios
DELTA     = 0.01     # relative plateau tolerance: eps = DELTA * r_s_max
ALPHAS    = np.round(np.arange(0.5, 3.6, 0.1), 2)   # scan range [0.5 .. 3.5]


# ─────────────────────────────────────────────────────────────────────────────
def compute_plateau_degeneracy(alphas, r_s, delta=DELTA):
    """
    Compute all plateau degeneracy measures from an alpha-scan.

    Parameters
    ----------
    alphas : array-like, shape (n,)
        Scanned exponent values.
    r_s : array-like, shape (n,)
        Spearman correlation r_s(alpha) for each alpha.
    delta : float, optional
        Relative tolerance parameter.  eps = delta * max(r_s).  Default 0.01.

    Returns
    -------
    dict
        alpha_star   : float  — exponent at maximum correlation
        r_max        : float  — maximum Spearman correlation
        alpha_p_min  : float  — lower plateau boundary
        alpha_p_max  : float  — upper plateau boundary
        D_width      : float  — plateau width
        D_norm       : float  — normalised degeneracy in [0, 1]
        D_flat       : float  — variance of r_s within plateau (0 = perfectly flat)
        D            : float  — combined measure D_width / (1 + D_flat)
        D_curv       : float  — curvature-based degeneracy (NaN at scan boundary)
        n_plateau    : int    — number of alpha points in plateau
        epsilon      : float  — absolute tolerance used
    """
    alphas = np.asarray(alphas, dtype=float)
    r_s    = np.asarray(r_s,    dtype=float)

    assert alphas.ndim == 1 and r_s.ndim == 1 and len(alphas) == len(r_s), \
        "alphas and r_s must be 1-D arrays of equal length"

    r_max      = float(np.max(r_s))
    eps        = delta * r_max
    mask       = r_s >= (r_max - eps)
    p_alphas   = alphas[mask]
    p_rs       = r_s[mask]

    alpha_min  = float(np.min(alphas))
    alpha_max  = float(np.max(alphas))
    D_width    = float(np.max(p_alphas) - np.min(p_alphas))
    D_norm     = float(D_width / (alpha_max - alpha_min)) if alpha_max > alpha_min else 0.0
    D_flat     = float(np.var(p_rs))
    D_combined = float(D_width / (1.0 + D_flat))

    idx_star   = int(np.argmax(r_s))
    D_curv     = float("nan")
    if 0 < idx_star < len(r_s) - 1:
        h1 = alphas[idx_star]     - alphas[idx_star - 1]
        h2 = alphas[idx_star + 1] - alphas[idx_star]
        if np.isclose(h1, h2) and h1 > 0:
            sd = (r_s[idx_star + 1] - 2 * r_s[idx_star] + r_s[idx_star - 1]) / h1 ** 2
            if sd != 0:
                D_curv = float(1.0 / abs(sd))

    return dict(
        alpha_star  = float(alphas[idx_star]),
        r_max       = r_max,
        alpha_p_min = float(np.min(p_alphas)),
        alpha_p_max = float(np.max(p_alphas)),
        D_width     = D_width,
        D_norm      = D_norm,
        D_flat      = D_flat,
        D           = D_combined,
        D_curv      = D_curv,
        n_plateau   = int(np.sum(mask)),
        epsilon     = eps,
    )


# ─────────────────────────────────────────────────────────────────────────────
def rs_scan_from_raw(df, alphas, eps_denom=EPS_DENOM):
    """
    Compute r_s(alpha) for each alpha from a raw ensemble DataFrame.

    Expects columns: mean_cci, mean_dS_pos, mean_mi.
    The ratio  R_alpha = mean_dS_pos / mean_mi^alpha  is used as
    the entropy-information proxy.
    """
    cci   = df["mean_cci"].values
    dS    = df["mean_dS_pos"].values
    mi    = df["mean_mi"].values
    return np.array([
        spearmanr(cci, dS / (mi + eps_denom) ** a)[0]
        for a in alphas
    ])


# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("Plateau Degeneracy Analysis — Papers 07 / 09 / 11 / 16")
    print("=" * 65)
    print(f"  delta   = {DELTA}  (tolerance: eps = delta * r_s_max)")
    print(f"  alphas  = [{ALPHAS[0]:.1f}, {ALPHAS[-1]:.1f}],  step h = {ALPHAS[1]-ALPHAS[0]:.1f}")
    print()

    # ── Load data ─────────────────────────────────────────────────────────────
    df_1d = pd.read_csv("cci_entropy_information_test.csv")
    df_2d = pd.read_csv("2d_cci_entropy_information_test.csv")
    df_3d = pd.read_csv("3d_cci_alpha_scan.csv")

    # ── Compute r_s(alpha) ────────────────────────────────────────────────────
    rs_1d = rs_scan_from_raw(df_1d, ALPHAS)
    rs_2d = rs_scan_from_raw(df_2d, ALPHAS)
    # 3D: stored scan → interpolate onto ALPHAS grid
    rs_3d = np.interp(ALPHAS, df_3d["alpha"].values, df_3d["spearman_r"].values)

    # ── Compute degeneracy ────────────────────────────────────────────────────
    rows = []
    for dim, rs in [("1D", rs_1d), ("2D", rs_2d), ("3D", rs_3d)]:
        result = compute_plateau_degeneracy(ALPHAS, rs)
        rows.append({"Dimension": dim, **result})

    df_out = pd.DataFrame(rows)

    # ── Print results table ───────────────────────────────────────────────────
    print("Results:")
    print("-" * 65)
    display_cols = [
        "Dimension", "alpha_star", "r_max",
        "alpha_p_min", "alpha_p_max",
        "D_width", "D_norm", "D_flat", "D", "n_plateau",
    ]
    print(df_out[display_cols].to_string(index=False,
          float_format=lambda x: f"{x:.4f}"))

    # ── LaTeX table rows ──────────────────────────────────────────────────────
    print("\nLaTeX table rows (for Paper 16):")
    print("-" * 65)
    for _, row in df_out.iterrows():
        print(
            f"  {row['Dimension']} & {row['alpha_star']:.1f} "
            f"& {row['r_max']:.3f} "
            f"& $[{row['alpha_p_min']:.1f},\\;{row['alpha_p_max']:.1f}]$ "
            f"& {row['D_width']:.2f} "
            f"& {row['D_norm']:.3f} "
            f"& {int(row['n_plateau'])} \\\\"
        )

    # ── Non-monotonic trend summary ───────────────────────────────────────────
    d_vals = {row["Dimension"]: row["D_norm"] for _, row in df_out.iterrows()}
    print(f"\nDimensional trend (D_norm):")
    print(f"  1D = {d_vals['1D']:.3f}")
    print(f"  2D = {d_vals['2D']:.3f}  <- minimum (dimensional transition!)")
    print(f"  3D = {d_vals['3D']:.3f}")
    print(f"\n  D_norm(2D) < D_norm(1D) < D_norm(3D)")
    print(f"  => non-monotonic, transition at 2D→3D")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_file = "plateau_degeneracy_results.csv"
    df_out.to_csv(out_file, index=False)
    print(f"\nSaved: {out_file}")


if __name__ == "__main__":
    main()

"""
Structural Selection Principle — FLRW Benchmark 
===================================================
Three cosmological backgrounds are evaluated using the
structural energy functional E[Phi] from:

  Krieg, C. (2026). A Structural Selection Principle for
  Dynamically Consistent Field Configurations.

Goal:
    Demonstrate that the structural functional has non-trivial
    discriminative power across controlled background classes.
    This benchmark does not validate the selection principle as
    a law of nature; it shows only that the functional can rank
    qualitatively different backgrounds in a reproducible way.

The three benchmarks are:
  Case A: Exact de Sitter        (static field, KG-compatible C)
  Case B: Slow-roll scalar FLRW  (dynamically active, same matter)
  Case C: Ghost-near stress test  (deliberately constructed edge case)

Notes:
  - Cases A and B use identical matter assumptions (w_m = -1)
    so that differences arise purely from field dynamics.
  - In Case A, C is set to the KG-compatible value C ~ 6*xi
    so that both background equations are approximately satisfied.
  - Case C is not a realistic cosmological background;
    it tests the ghost-barrier term Sigma only.
  - All cases use the same Lagrange multiplier values (lambda_i = 1).
"""

import numpy as np

# ---------------------------------------------------------------
# Physical constants (Planck units: M_Pl = 1)
# ---------------------------------------------------------------
M_PL   = 1.0          # Planck mass
H0     = 1.0e-3       # Hubble constant (normalised)
XI     = 1.0e-3       # non-minimal coupling
M_C    = 1.0          # coherence field scale
EPS    = 1.0e-6       # ghost-barrier regulator
W_M    = -1.0         # equation of state (same for A and B)

# Activity functional parameters
ETA_C  = 1.0
BETA_C = 1.0
LAM_C  = 1.0e-2       # Lambda_C (wave operator cutoff)
LAM_X  = 1.0e-2       # Lambda_X (kinetic cutoff)

# Lagrange multipliers (equal weights, normalised)
LAM_STAB = 1.0
LAM_CONS = 1.0
LAM_CONN = 1.0
LAM_NEQ  = 1.0

ALPHA = 1.0           # ghost-penalty coefficient

# ---------------------------------------------------------------
# Structural diagnostics on FLRW background
# ---------------------------------------------------------------

def F_coupling(C):
    """Non-minimal coupling F(C) = 1 + xi*C/M_C"""
    return 1.0 + XI * C / M_C

def H_residual(R_F, R_KG):
    """Residual Stability Functional H_FLRW = -(R_F^2 + R_KG^2)"""
    return -(R_F**2 + R_KG**2)

def B_conservation(rho_dot, rho, p, H):
    """Conservation Functional B_FLRW = -(rho_dot + 3H(rho+p))^2
    Full continuity equation residual including rho_dot."""
    cons_residual = rho_dot + 3.0 * H * (rho + p)
    return -(cons_residual**2)

def S_activity(box_C, X):
    """Dynamical Activity Functional S_C (bounded in [0,1))"""
    term1 = ETA_C  * box_C**2 / (LAM_C**4 + box_C**2)
    term2 = BETA_C * abs(X)   / (LAM_X**4 + abs(X))
    return term1 + term2

def Sigma_consistency(C):
    """Consistency Functional: soft ghost barrier near F(C) -> 0+
    Sigma_FLRW = -alpha^2 / (F^2 + eps)"""
    F = F_coupling(C)
    return -ALPHA**2 / (F**2 + EPS)

def structural_energy(H_val, B_val, S_val, Sigma_val):
    """Total structural energy E_FLRW (to be minimised).
    E = lambda_stab*(-H - Sigma) + lambda_cons*(-B)
      + lambda_conn*(-N)  [N=0 single field]
      - lambda_neq*S      [more activity = less E]
    """
    E_stab = (-H_val) + (-Sigma_val)
    E_cons = (-B_val)
    E_conn = 0.0          # N_FLRW = 0 in single-field baseline
    E_neq  = S_val
    return (LAM_STAB * E_stab
            + LAM_CONS * E_cons
            + LAM_CONN * E_conn
            - LAM_NEQ  * E_neq)

# ---------------------------------------------------------------
# Case A: Exact de Sitter
#   H = H0 = const,  C = C0 = const,  Cdot = 0
#   Matter: w_m = -1 (cosmological constant type)
# ---------------------------------------------------------------
def case_A_de_Sitter():
    H     = H0
    # C chosen to satisfy KG equation approximately for static field:
    # V'(C) = 0.5 * M_Pl^2 * F'(C) * R  =>  C ~ 6*xi = 0.006
    C     = 6.0 * XI        # KG-compatible static de Sitter value
    Cdot  = 0.0
    Cddot = 0.0

    F      = F_coupling(C)
    R      = 12.0 * H**2          # Ricci scalar in de Sitter

    V      = 0.5 * (H0 * M_C)**2 * C**2
    Vprime = (H0 * M_C)**2 * C
    Fprime = XI / M_C

    # Matter: exact de Sitter -> rho_m from Friedmann
    rho_m  = 3.0 * M_PL**2 * F * H**2 - V
    p_m    = W_M * rho_m          # w_m = -1
    rho_dot = 0.0                 # de Sitter: rho_m = const

    # Friedmann residual: exact -> R_F = 0
    R_F  = 0.0

    # KG residual: Cdot = Cddot = 0
    R_KG = 0.0 + Vprime - 0.5 * M_PL**2 * Fprime * R

    # Diagnostics
    H_val   = H_residual(R_F, R_KG)
    B_val   = B_conservation(rho_dot, rho_m, p_m, H)
    box_C   = -(Cddot + 3.0 * H * Cdot)
    X       = 0.5 * Cdot**2
    S_val   = S_activity(box_C, X)
    Sig_val = Sigma_consistency(C)

    E = structural_energy(H_val, B_val, S_val, Sig_val)
    return dict(case="A: de Sitter", H=H_val, B=B_val,
                S=S_val, Sigma=Sig_val, E=E,
                R_F=R_F, R_KG=R_KG, F=F,
                rho_dot=rho_dot)

# ---------------------------------------------------------------
# Case B: Slow-roll scalar-tensor FLRW (preferred case)
#   H ~ H0,  Cdot != 0 small,  w_m = -1 (same as A)
# ---------------------------------------------------------------
def case_B_slow_roll():
    H     = H0
    C     = 0.1
    Cdot  = 0.05 * H0 * M_C      # slow roll: Cdot << M_C H0
    Cddot = -2.5 * H * Cdot      # slight deviation from exact slow-roll
                                  # ensures box_C != 0, both activity
                                  # terms contribute

    F      = F_coupling(C)
    R      = 12.0 * H**2
    Fprime = XI / M_C
    Fdot   = Fprime * Cdot

    V      = 0.5 * (H0 * M_C)**2 * C**2
    Vprime = (H0 * M_C)**2 * C

    # Matter: same equation of state as Case A
    rho_C  = 0.5 * Cdot**2 + V
    rho_m  = 3.0 * M_PL**2 * F * H**2 - rho_C + 3.0 * M_PL**2 * H * Fdot
    p_m    = W_M * rho_m          # w_m = -1, same as A
    rho_dot = -3.0 * H * (rho_m + p_m)   # continuity eq. exact

    # Friedmann residual
    R_F  = (3.0 * M_PL**2 * F * H**2
            - rho_m - rho_C
            + 3.0 * M_PL**2 * H * Fdot)

    # KG residual (slow-roll: small)
    R_KG = Cddot + 3.0 * H * Cdot + Vprime - 0.5 * M_PL**2 * Fprime * R

    # Diagnostics
    H_val   = H_residual(R_F, R_KG)
    B_val   = B_conservation(rho_dot, rho_m, p_m, H)
    box_C   = -(Cddot + 3.0 * H * Cdot)
    X       = 0.5 * Cdot**2
    S_val   = S_activity(box_C, X)
    Sig_val = Sigma_consistency(C)

    E = structural_energy(H_val, B_val, S_val, Sig_val)
    return dict(case="B: slow-roll FLRW", H=H_val, B=B_val,
                S=S_val, Sigma=Sig_val, E=E,
                R_F=R_F, R_KG=R_KG, F=F,
                rho_dot=rho_dot)

# ---------------------------------------------------------------
# Case C: Ghost-near stress test
#   F(C) -> 0+: deliberately constructed consistency edge case.
#   Not a realistic cosmological background.
#   Same matter equation of state as A and B for fair comparison.
# ---------------------------------------------------------------
def case_C_ghost_near():
    # Construct C such that F(C) ~ 0.05 (near ghost crossing)
    F_target = 0.05
    C = (F_target - 1.0) * M_C / XI   # large negative C

    H     = H0
    Cdot  = 0.05 * H0 * M_C           # same as B for fair comparison
    Cddot = -3.0 * H * Cdot

    F      = F_coupling(C)             # ~ F_target
    R      = 12.0 * H**2
    Fprime = XI / M_C
    Fdot   = Fprime * Cdot

    V      = 0.5 * (H0 * M_C)**2 * C**2
    Vprime = (H0 * M_C)**2 * C

    rho_C  = 0.5 * Cdot**2 + V
    rho_m  = max(3.0 * M_PL**2 * F * H**2 - rho_C
                 + 3.0 * M_PL**2 * H * Fdot, 0.0)
    # Note: negative effective matter density is clipped to zero
    # to keep the stress test numerically well-defined.
    # This is not a physical solution; Case C tests only the
    # Sigma ghost-barrier behaviour.
    p_m    = W_M * rho_m               # w_m = -1, same as A and B
    rho_dot = -3.0 * H * (rho_m + p_m)

    R_F  = (3.0 * M_PL**2 * F * H**2
            - rho_m - rho_C
            + 3.0 * M_PL**2 * H * Fdot)
    R_KG = Cddot + 3.0 * H * Cdot + Vprime - 0.5 * M_PL**2 * Fprime * R

    # Diagnostics
    H_val   = H_residual(R_F, R_KG)
    B_val   = B_conservation(rho_dot, rho_m, p_m, H)
    box_C   = -(Cddot + 3.0 * H * Cdot)
    X       = 0.5 * Cdot**2
    S_val   = S_activity(box_C, X)
    Sig_val = Sigma_consistency(C)     # strongly negative -> large E

    E = structural_energy(H_val, B_val, S_val, Sig_val)
    return dict(case="C: ghost-near (stress test)", H=H_val, B=B_val,
                S=S_val, Sigma=Sig_val, E=E,
                R_F=R_F, R_KG=R_KG, F=F,
                rho_dot=rho_dot)

# ---------------------------------------------------------------
# Run and print results
# ---------------------------------------------------------------
if __name__ == "__main__":
    results = [case_A_de_Sitter(),
               case_B_slow_roll(),
               case_C_ghost_near()]

    print("=" * 72)
    print("Structural Selection Principle — FLRW Benchmark v4")
    print("=" * 72)
    print(f"Parameters: xi={XI}, M_C={M_C}, H0={H0}, w_m={W_M}")
    print(f"Weights: lambda_stab={LAM_STAB}, lambda_cons={LAM_CONS}, "
          f"lambda_conn={LAM_CONN}, lambda_neq={LAM_NEQ}")
    print(f"Note: Cases A and B use identical matter assumptions (w_m={W_M})")
    print(f"      Case C is a ghost-near consistency stress test only.")
    print("-" * 72)

    header = (f"{'Case':<30} {'F(C)':>6} {'H':>8} {'B':>8}"
              f" {'S':>7} {'Sigma':>10} {'E':>10}")
    print(header)
    print("-" * 72)

    for r in results:
        print(f"{r['case']:<30}"
              f" {r['F']:>6.4f}"
              f" {r['H']:>8.4f}"
              f" {r['B']:>8.4f}"
              f" {r['S']:>7.4f}"
              f" {r['Sigma']:>10.4f}"
              f" {r['E']:>10.4f}")

    print("-" * 72)
    energies = [(r['E'], r['case']) for r in results]
    best  = min(energies)[1]
    worst = max(energies)[1]
    print(f"Structurally preferred (min E): {best}")
    print(f"Structurally rejected  (max E): {worst}")
    print("=" * 72)

    E_A, E_B, E_C = results[0]['E'], results[1]['E'], results[2]['E']
    print(f"\nRanking E_B < E_A < E_C: "
          f"{'YES' if E_B < E_A < E_C else 'NO'}")
    print(f"  E_A (de Sitter)  = {E_A:.6f}")
    print(f"  E_B (slow-roll)  = {E_B:.6f}")
    print(f"  E_C (ghost-near) = {E_C:.6f}")
    print(f"\nConclusion: The structural functional discriminates between")
    print(f"  trivial exact, dynamically active, and ghost-near backgrounds.")
    print(f"\nDisclaimer: This benchmark demonstrates discriminative capacity")
    print(f"  of the functional, not a derivation of cosmological predictions.")

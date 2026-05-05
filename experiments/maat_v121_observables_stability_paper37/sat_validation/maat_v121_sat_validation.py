#!/usr/bin/env python3
# maat_v121_sat_validation.py

import argparse
import json
import os
import numpy as np
import pandas as pd

try:
    from scipy.stats import pearsonr, spearmanr
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False

try:
    from sklearn.linear_model import Ridge, LinearRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    SKLEARN_OK = True
except Exception:
    SKLEARN_OK = False


EPS = 1e-12


def maat_v121_closure(H, B, S, V, eps=EPS):
    H = np.clip(np.asarray(H, dtype=float), eps, 1.0)
    B = np.clip(np.asarray(B, dtype=float), eps, 1.0)
    S = np.clip(np.asarray(S, dtype=float), eps, 1.0)
    V = np.clip(np.asarray(V, dtype=float), eps, 1.0)

    R_imp = np.power(H * B * V, 1.0 / 3.0)
    geom4 = np.power(H * B * S * V, 1.0 / 4.0)
    R_em = np.minimum(R_imp, geom4)

    CCI_v121 = S / (H + B + V + R_imp + R_em + eps)
    F_v121 = -(np.log(eps + H) + np.log(eps + B) + np.log(eps + S)
               + np.log(eps + V) + np.log(eps + R_em))

    return R_imp, R_em, CCI_v121, F_v121


def find_runtime_column(df):
    candidates = [
        "runtime", "runtime_s", "solve_time", "time", "seconds",
        "wall_time", "wall_clock_time", "elapsed", "elapsed_s",
        "log_runtime", "log10_runtime", "median_runtime",
        "median_log_runtime", "log10_time"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return None


def add_fields_from_available_columns(df):
    # Case 1: already field supports
    if all(c in df.columns for c in ["H", "B", "S", "V"]):
        df["H"] = df["H"].astype(float)
        df["B"] = df["B"].astype(float)
        df["S"] = df["S"].astype(float)
        df["V"] = df["V"].astype(float)
        if "R" in df.columns:
            df["R"] = df["R"].astype(float)
        return df, "field_columns_HBSV"

    # Case 2: defect style dH
    if all(c in df.columns for c in ["dH", "dB", "dS", "dV"]):
        df["H"] = 1.0 / (1.0 + df["dH"].astype(float))
        df["B"] = 1.0 / (1.0 + df["dB"].astype(float))
        df["S"] = 1.0 / (1.0 + df["dS"].astype(float))
        df["V"] = 1.0 / (1.0 + df["dV"].astype(float))
        if "dR" in df.columns:
            df["R"] = 1.0 / (1.0 + df["dR"].astype(float))
        return df, "defect_columns_dH_dB_dS_dV"

    # Case 3: defect style d_H
    if all(c in df.columns for c in ["d_H", "d_B", "d_S", "d_V"]):
        df["H"] = 1.0 / (1.0 + df["d_H"].astype(float))
        df["B"] = 1.0 / (1.0 + df["d_B"].astype(float))
        df["S"] = 1.0 / (1.0 + df["d_S"].astype(float))
        df["V"] = 1.0 / (1.0 + df["d_V"].astype(float))
        if "d_R" in df.columns:
            df["R"] = 1.0 / (1.0 + df["d_R"].astype(float))
        return df, "defect_columns_d_H_d_B_d_S_d_V"

    raise ValueError(
        "Missing columns. Need H,B,S,V or dH,dB,dS,dV or d_H,d_B,d_S,d_V.\n"
        f"Available columns: {list(df.columns)}"
    )


def print_group_summary(df, group_col):
    if group_col not in df.columns:
        return None

    cols = ["H", "B", "S", "V", "R_imp_v121", "R_em_v121", "CCI_v121", "F_v121"]
    summary = df.groupby(group_col)[cols].agg(["mean", "std", "min", "max"])

    print(f"\n--- Group summary by {group_col} ---")
    print(summary.round(6).to_string())

    out_name = f"maat_v121_group_summary_by_{group_col}.csv"
    summary.to_csv(out_name)
    return out_name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Input CSV file.")
    parser.add_argument("--out", default="maat_v121_sat_validation_results.csv")
    args = parser.parse_args()

    if not os.path.exists(args.csv):
        print(f"\n❌ File not found: {args.csv}")
        print("\nAvailable CSV files here:")
        for f in os.listdir("."):
            if f.endswith(".csv"):
                print(" -", f)
        raise SystemExit(1)

    df = pd.read_csv(args.csv)
    df, input_mode = add_fields_from_available_columns(df)

    H = df["H"].to_numpy(float)
    B = df["B"].to_numpy(float)
    S = df["S"].to_numpy(float)
    V = df["V"].to_numpy(float)

    R_imp, R_em, CCI_v121, F_v121 = maat_v121_closure(H, B, S, V)

    df["R_imp_v121"] = R_imp
    df["R_em_v121"] = R_em
    df["CCI_v121"] = CCI_v121
    df["F_v121"] = F_v121

    if "R" in df.columns:
        R_old = np.clip(df["R"].to_numpy(float), EPS, 1.0)
        legacy_stability = np.minimum(
            R_old,
            np.power(np.clip(H * B * S * V, EPS, 1.0), 0.25)
        )
        df["legacy_stability"] = legacy_stability
        df["delta_R_em_minus_legacy"] = R_em - legacy_stability

    print("\n=== MAAT v1.2.1 SAT / Defect Validation ===")
    print(f"Input file : {args.csv}")
    print(f"Rows       : {len(df)}")
    print(f"Input mode : {input_mode}")

    runtime_col = find_runtime_column(df)
    summary = {
        "input_file": args.csv,
        "rows": int(len(df)),
        "input_mode": input_mode,
        "runtime_column": runtime_col,
    }

    print("\n--- Overall v1.2.1 metrics ---")
    metric_cols = ["H", "B", "S", "V", "R_imp_v121", "R_em_v121", "CCI_v121", "F_v121"]
    print(df[metric_cols].describe().round(6).to_string())

    saved_group_files = []
    for group_col in ["label", "source"]:
        g = print_group_summary(df, group_col)
        if g:
            saved_group_files.append(g)

    if runtime_col is None:
        print("\n⚠️ No runtime column found.")
        print("This file looks like a defect/calibration dataset, not a runtime dataset.")
        print("The script still computed v1.2.1 closure metrics.")
        print("For runtime validation, use a CSV from results/ that contains runtime/time/seconds.")
    else:
        runtime_raw = df[runtime_col].astype(float).to_numpy()
        if runtime_col.startswith("log") or "log" in runtime_col:
            log_runtime = runtime_raw
        else:
            log_runtime = np.log10(np.maximum(runtime_raw, EPS))

        df["log10_runtime_used"] = log_runtime

        print(f"\nRuntime column found: {runtime_col}")

        if SCIPY_OK:
            print("\n--- Correlation with log10(runtime) ---")
            rows = []
            for name in ["H", "B", "S", "V", "R_imp_v121", "R_em_v121", "CCI_v121", "F_v121"]:
                values = df[name].to_numpy(float)
                pr, pp = pearsonr(values, log_runtime)
                sr, sp = spearmanr(values, log_runtime)
                rows.append({
                    "metric": name,
                    "pearson_r": pr,
                    "pearson_p": pp,
                    "spearman_rho": sr,
                    "spearman_p": sp,
                })
            corr_df = pd.DataFrame(rows)
            print(corr_df.round(6).to_string(index=False))
            corr_df.to_csv("maat_v121_sat_correlations.csv", index=False)

        if SKLEARN_OK and len(df) >= 10:
            X_base_cols = ["H", "B", "S", "V"]
            X_v121_cols = ["H", "B", "S", "V", "R_imp_v121", "R_em_v121", "CCI_v121"]

            if "alpha" in df.columns:
                X_base_cols = ["alpha"] + X_base_cols
                X_v121_cols = ["alpha"] + X_v121_cols

            X_base = df[X_base_cols].to_numpy(float)
            X_v121 = df[X_v121_cols].to_numpy(float)

            cv = min(5, len(df))
            model_base = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
            model_v121 = make_pipeline(StandardScaler(), Ridge(alpha=1.0))

            cv_base = cross_val_score(model_base, X_base, log_runtime, cv=cv, scoring="r2")
            cv_v121 = cross_val_score(model_v121, X_v121, log_runtime, cv=cv, scoring="r2")

            print("\n--- Ridge CV R² ---")
            print(f"Base fields    : mean={cv_base.mean():.6f}, std={cv_base.std():.6f}")
            print(f"v1.2.1 closure : mean={cv_v121.mean():.6f}, std={cv_v121.std():.6f}")

            summary["base_cv_r2_mean"] = float(cv_base.mean())
            summary["v121_cv_r2_mean"] = float(cv_v121.mean())
            summary["delta_cv_r2"] = float(cv_v121.mean() - cv_base.mean())

            lin = make_pipeline(StandardScaler(), LinearRegression())
            lin.fit(X_v121, log_runtime)
            coefs = lin.named_steps["linearregression"].coef_

            coef_df = pd.DataFrame({
                "feature": X_v121_cols,
                "standardized_coef": coefs,
            }).sort_values("standardized_coef", key=lambda s: np.abs(s), ascending=False)

            print("\n--- Standardized coefficients ---")
            print(coef_df.round(6).to_string(index=False))
            coef_df.to_csv("maat_v121_sat_coefficients.csv", index=False)

    df.to_csv(args.out, index=False)

    with open("maat_v121_sat_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\nSaved:")
    print(f"- {args.out}")
    print("- maat_v121_sat_summary.json")
    for f in saved_group_files:
        print(f"- {f}")


if __name__ == "__main__":
    main()
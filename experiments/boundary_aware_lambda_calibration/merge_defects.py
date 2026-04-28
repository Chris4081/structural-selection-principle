# merge_defects

import pandas as pd

sat = pd.read_csv("maat_defects_sat.csv")
core = pd.read_csv("maat_defects_core.csv")
boundary = pd.read_csv("maat_defects_core_boundary.csv")

cols = ["d_H", "d_B", "d_S", "d_V", "d_R", "label", "source"]

df = pd.concat(
    [
        sat[cols],
        core[cols],
        boundary[cols],
    ],
    ignore_index=True,
)

df.to_csv("maat_defects_fused.csv", index=False)

print("Wrote maat_defects_fused.csv")
print("\n--- Source counts ---")
print(df["source"].value_counts())

print("\n--- Label counts ---")
print(df["label"].value_counts())

print("\n--- Defect summary ---")
print(df[["d_H", "d_B", "d_S", "d_V", "d_R"]].describe())
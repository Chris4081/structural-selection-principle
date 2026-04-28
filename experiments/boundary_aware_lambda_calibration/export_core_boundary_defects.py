import numpy as np
import pandas as pd
from maat_core import MaatCore, Field, Constraint

rng = np.random.default_rng(148)
rows = []

SECTORS = ["H", "B", "S", "V", "R"]

def state_fn(arr):
    return {k: float(v) for k, v in zip(SECTORS, arr)}

def make_field(name, target, weight=1.0):
    return Field(
        name=name,
        func=lambda x, n=name, t=target: (x[n] - t) ** 2,
        weight=float(weight),
    )

def run_case(label, i):
    if label == "safe":
        thresholds = {"H": 0.35, "B": 0.35, "V": 0.20, "R": 0.45}
        targets = {"H": 0.80, "B": 0.80, "S": 0.70, "V": 0.75, "R": 0.85}
        safety_lambda = 1e6
    elif label == "near_boundary":
        thresholds = {"H": 0.72, "B": 0.72, "V": 0.68, "R": 0.78}
        targets = {"H": 0.76, "B": 0.76, "S": 0.65, "V": 0.71, "R": 0.80}
        safety_lambda = 1e5
    else:  # violated stress test
        thresholds = {"H": 0.92, "B": 0.92, "V": 0.92, "R": 0.95}
        targets = {"H": 0.55, "B": 0.55, "S": 0.40, "V": 0.55, "R": 0.60}
        safety_lambda = 10.0

    targets = {k: float(np.clip(v + rng.normal(0, 0.03), 0.0, 1.0)) for k, v in targets.items()}

    fields = [
        make_field("H", targets["H"], rng.uniform(0.8, 1.8)),
        make_field("B", targets["B"], rng.uniform(0.8, 1.8)),
        make_field("S", targets["S"], rng.uniform(0.8, 1.8)),
        make_field("V", targets["V"], rng.uniform(0.8, 1.8)),
        make_field("R", targets["R"], rng.uniform(0.8, 1.8)),
    ]

    constraints = [
        Constraint("H_min", lambda x, t=thresholds["H"]: x["H"] - t),
        Constraint("B_min", lambda x, t=thresholds["B"]: x["B"] - t),
        Constraint("V_min", lambda x, t=thresholds["V"]: x["V"] - t),
        Constraint("R_min", lambda x, t=thresholds["R"]: x["R"] - t),
        Constraint("upper_H", lambda x: 1.0 - x["H"]),
        Constraint("upper_B", lambda x: 1.0 - x["B"]),
        Constraint("upper_S", lambda x: 1.0 - x["S"]),
        Constraint("upper_V", lambda x: 1.0 - x["V"]),
        Constraint("upper_R", lambda x: 1.0 - x["R"]),
        Constraint("lower_H", lambda x: x["H"]),
        Constraint("lower_B", lambda x: x["B"]),
        Constraint("lower_S", lambda x: x["S"]),
        Constraint("lower_V", lambda x: x["V"]),
        Constraint("lower_R", lambda x: x["R"]),
    ]

    core = MaatCore(fields=fields, constraints=constraints, safety_lambda=safety_lambda)

    x0 = rng.uniform(0.0, 1.0, size=5)

    result = core.seek(
        state_fn,
        x0,
        bounds=((0.0, 1.0),) * 5,
        maxiter=1000,
        method="L-BFGS-B",
        seed=1000 + i,
        scalar_compat=False,
    )

    x = state_fn(result.x)
    margins = np.array([float(c.func(x)) for c in constraints])
    min_margin = float(margins.min())
    violation = max(0.0, -min_margin)

    # Respect defect: not just violation, also boundary stress
    boundary_stress = float(np.exp(-max(min_margin, 0.0) / 0.05))
    d_R = float(np.clip(violation + 0.25 * boundary_stress, 0.0, 1.0))

    rows.append({
        "d_H": 1.0 - x["H"],
        "d_B": 1.0 - x["B"],
        "d_S": 1.0 - x["S"],
        "d_V": 1.0 - x["V"],
        "d_R": d_R,
        "label": label,
        "source": "maat_core_boundary",
        "min_margin": min_margin,
        "violation": violation,
        "boundary_stress": boundary_stress,
        "safety_lambda": safety_lambda,
    })

for label in ["safe", "near_boundary", "violated"]:
    for i in range(300):
        run_case(label, i)

df = pd.DataFrame(rows)
df.to_csv("../maat_defects_core_boundary.csv", index=False)

print("Wrote ../maat_defects_core_boundary.csv")
print(df["label"].value_counts())
print(df.groupby("label")[["d_H", "d_B", "d_S", "d_V", "d_R", "min_margin", "violation"]].mean())
print(df.describe())

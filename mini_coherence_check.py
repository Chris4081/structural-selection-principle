
import numpy as np
import matplotlib.pyplot as plt

# 0D coherence dynamics
a = 1.0
b = 1.0
dt = 0.01
T = 12.0
times = np.arange(0, T + dt, dt)

def simulate_0d(C0):
    C = C0
    vals, energies = [], []
    for _ in times:
        vals.append(C)
        energies.append(-0.5 * a * C**2 + 0.25 * b * C**4)
        C += dt * (a * C - b * C**3)
    return np.array(vals), np.array(energies)

# 1D coherence field
D = 0.8
L = 128
dx = 1.0
dt1 = 0.05
steps = 600
rng = np.random.default_rng(42)
C = 0.05 * rng.normal(size=L)

def laplacian_periodic(u):
    return (np.roll(u, -1) - 2*u + np.roll(u, 1)) / dx**2

def free_energy_1d(u):
    grad = (np.roll(u, -1) - u) / dx
    density = 0.5 * D * grad**2 - 0.5 * a * u**2 + 0.25 * b * u**4
    return density.sum() * dx

energies = []
for _ in range(steps):
    energies.append(free_energy_1d(C))
    C = C + dt1 * (D * laplacian_periodic(C) + a * C - b * C**3)

print("0D stable fixed points: ±", np.sqrt(a/b))
print("1D final field mean:", C.mean())
print("1D final field std:", C.std())
print("1D initial energy > final energy:", energies[0], ">", energies[-1])

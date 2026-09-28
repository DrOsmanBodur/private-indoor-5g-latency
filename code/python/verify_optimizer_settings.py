"""Checks the optimizer settings of Section 4.3.2 without MATLAB.

A bounded Levenberg-Marquardt fit is started from the initial vector psi_0
reported in the article and, for the multi-start analysis, from each reported
initial value of C:

    psi_0 = [8.64, 0.0305, 486.89, 611.39, 0, ..., 0]
    0 <= L0 <= 50, 0 <= beta <= 5, 0 <= theta <= 5000,
    max(b) + 1 <= C <= 5000, -50 <= gamma_q <= 50
    C_0 in {max(b) + 1, 530, 611.39, 800, 1200, 2000, 5000} Mbps

Every start should reach the published estimates.
"""
import csv
import os

import numpy as np

# Paths are resolved relative to this script so the package is self-contained.
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(_HERE, "..", "..", "data"))

lat, bw = {}, {}
for r in csv.DictReader(open(os.path.join(DATA, "latency_210.csv"), encoding="utf-8-sig")):
    lat[(int(r["qIdx"]), int(float(r["B_Mbps"])), int(r["trial"]))] = float(r["latency_ms"])
for r in csv.DictReader(open(os.path.join(DATA, "latency_trials_raw_210.csv"),
                             encoding="utf-8")):
    bw[(int(r["qIdx"]), int(r["target_Mbps"]), int(r["trial"]))] = float(r["achieved_Mbps"])
keys = sorted(lat)
q = np.array([k[0] for k in keys]); L = np.array([lat[k] for k in keys])
b = np.array([bw[k] for k in keys])
Q = np.zeros((len(keys), 6))
for j in range(2, 8):
    Q[:, j - 2] = (q == j)


def model(p):
    L0, be, th, C = p[:4]
    return L0 + be * b + th / (C - b) + Q @ p[4:]


def jac(p):
    L0, be, th, C = p[:4]
    d = C - b
    J = np.zeros((len(b), 10))
    J[:, 0] = 1.0
    J[:, 1] = b
    J[:, 2] = 1.0 / d
    J[:, 3] = -th / d ** 2
    J[:, 4:] = Q
    return J


def fit(p0, lb, ub, max_iter=20000):
    """Bounded (projected) Levenberg-Marquardt; returns parameters and iterations."""
    p = np.clip(p0.copy(), lb, ub)
    r = L - model(p); s = float(r @ r); lam = 1e-3
    for it in range(1, max_iter + 1):
        J = jac(p); A = J.T @ J
        pn = np.clip(p + np.linalg.solve(A + lam * np.diag(np.diag(A)), J.T @ r), lb, ub)
        rn = L - model(pn); sn = float(rn @ rn)
        if sn < s:
            done = (s - sn) <= 1e-15 * s
            p, r, s, lam = pn, rn, sn, max(lam / 3, 1e-15)
            if done: break
        else:
            lam *= 4
            if lam > 1e14: break
    return p, it


bmax = b.max()
lb = np.array([0.0, 0.0, 0.0, bmax + 1.0] + [-50.0] * 6)
ub = np.array([50.0, 5.0, 5000.0, 5000.0] + [50.0] * 6)
psi0 = np.array([8.64, 0.0305, 486.89, 611.39] + [0.0] * 6)

print("Initial vector psi_0:", " ".join(f"{v:g}" for v in psi0))
print("Lower bounds        :", " ".join(f"{v:g}" for v in lb))
print("Upper bounds        :", " ".join(f"{v:g}" for v in ub))
print(f"  (lower bound of C = max(b) + 1 = {bmax + 1:.1f} Mbps)")
print()

p, it = fit(psi0, lb, ub)
r = L - model(p)
mae = np.mean(np.abs(r)); rmse = np.sqrt(np.mean(r ** 2))
r2 = 1 - np.sum(r ** 2) / np.sum((L - L.mean()) ** 2)
ref = [8.639, 0.0305, 486.885, 611.395, -0.470, -1.449, -1.306, 0.782, -0.476, -1.681]
nm = ["L0", "beta", "theta", "C", "g2", "g3", "g4", "g5", "g6", "g7"]
print(f"Fit from psi_0: {it} iterations")
print(f"  {'param':6s} {'fit':>12s} {'article':>12s} {'diff':>10s}")
for n, v, rr in zip(nm, p, ref):
    print(f"  {n:6s} {v:12.4f} {rr:12.4f} {abs(v-rr):10.4f}")
print(f"  MAE={mae:.3f} (article 2.52)  RMSE={rmse:.3f} (3.13)  R2={r2:.4f} (0.842)")
print(f"  parameters at a bound: {int(np.sum(np.isclose(p, lb) | np.isclose(p, ub)))}")

print()
print("Multi-start analysis over the initial value of C:")
print(f"  {'C0':>8s} {'iter':>6s} {'C':>10s} {'theta':>10s} {'SSE':>12s}")
for c0 in (bmax + 1, 530, 611.39, 800, 1200, 2000, 5000):
    start = psi0.copy(); start[3] = c0
    pc, itc = fit(start, lb, ub)
    sse = float(np.sum((L - model(pc)) ** 2))
    print(f"  {c0:8.2f} {itc:6d} {pc[3]:10.3f} {pc[2]:10.3f} {sse:12.6f}")

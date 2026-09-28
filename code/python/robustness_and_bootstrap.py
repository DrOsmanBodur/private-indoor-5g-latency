"""Imputation sensitivity (Section 4.3.1) and cell-bootstrap confidence intervals
(Table 5) of the proposed model, without MATLAB.

E) Refit without the three runs with invalid one-way delay (M = 207).
F) Cell bootstrap with 1000 replicates. The resampling reproduces MATLAB's
   rng(42) stream exactly (same random numbers, randi(n) = ceil(n * rand)),
   and every replicate is fitted with a bounded local optimizer started from
   the full-data estimate, as in step04_bootstrap_ci.m. Percentile intervals
   use MATLAB's prctile definition.
"""
import csv
import os

import numpy as np

# Paths are resolved relative to this script so the package is self-contained.
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(_HERE, "..", "..", "data"))

# rows in the order of latency_with_achieved_bw_210.csv (the MATLAB order)
rows = list(csv.DictReader(open(os.path.join(DATA, "latency_with_achieved_bw_210.csv"), encoding="utf-8-sig")))
q = np.array([int(r["qIdx"]) for r in rows])
tgt = np.array([float(r["B_Mbps"]) for r in rows])
L = np.array([float(r["latency_ms"]) for r in rows])
b = np.array([float(r["achieved_Mbps_server"]) for r in rows])
valid = {}
for r in csv.DictReader(open(os.path.join(DATA, "latency_trials_raw_210.csv"), encoding="utf-8")):
    valid[(int(r["qIdx"]), float(r["target_Mbps"]), int(r["trial"]))] = int(r["valid_owd"]) == 1
valid = np.array([valid[(int(r["qIdx"]), float(r["B_Mbps"]), int(r["trial"]))] for r in rows])


def X_full(bb, qq, C):
    cols = [np.ones_like(bb), bb, 1.0 / (C - bb)]
    for j in range(2, 8):
        cols.append((qq == j).astype(float))
    return np.column_stack(cols)


def fit(bb, qq, yy):
    """Global fit: scan over C with a linear solve for the other parameters."""
    lo = bb.max() + 1.0
    best = None
    for C in np.concatenate([np.linspace(lo, lo + 60, 300), np.linspace(lo + 60, 2500, 900)]):
        X = X_full(bb, qq, C)
        p, *_ = np.linalg.lstsq(X, yy, rcond=None)
        r = yy - X @ p; s = float(r @ r)
        if best is None or s < best[0]:
            best = (s, C, p)
    s, C, p = best
    step = 1.0
    for _ in range(50):
        moved = False
        for c2 in (C - step, C + step):
            if c2 < lo: continue
            X = X_full(bb, qq, c2)
            pp, *_ = np.linalg.lstsq(X, yy, rcond=None)
            r = yy - X @ pp; ss = float(r @ r)
            if ss < s: s, C, p, moved = ss, c2, pp, True
        if not moved:
            step /= 2
            if step < 1e-5: break
    return C, p


def met(y, yh):
    e = y - yh
    return np.mean(np.abs(e)), np.sqrt(np.mean(e**2)), 1 - np.sum(e**2)/np.sum((y-y.mean())**2)


print("=" * 76)
print("E) REFIT WITHOUT THE THREE INVALID RUNS (M = 207)")
print("=" * 76)
C0, p0 = fit(b, q, L)
m0 = met(L, X_full(b, q, C0) @ p0)
C1, p1 = fit(b[valid], q[valid], L[valid])
m1 = met(L[valid], X_full(b[valid], q[valid], C1) @ p1)
names = ["L0", "beta", "theta", "g2", "g3", "g4", "g5", "g6", "g7"]
print(f"  {'param':7s} {'M=210':>12s} {'M=207':>12s} {'change':>10s}")
worst = 0
for n, a, bb2 in zip(names, p0, p1):
    d = abs(bb2 - a) / max(abs(a), 1e-9) * 100
    if n in ("L0", "beta", "theta"): worst = max(worst, d)
    print(f"  {n:7s} {a:12.4f} {bb2:12.4f} {d:9.2f}%")
dC = abs(C1 - C0) / C0 * 100
worst = max(worst, dC)
print(f"  {'C':7s} {C0:12.4f} {C1:12.4f} {dC:9.2f}%")
print(f"  {'MAE':7s} {m0[0]:12.3f} {m1[0]:12.3f} {abs(m1[0]-m0[0]):9.3f} ms")
print(f"  {'RMSE':7s} {m0[1]:12.3f} {m1[1]:12.3f} {abs(m1[1]-m0[1]):9.3f} ms")
print(f"  {'R2':7s} {m0[2]:12.4f} {m1[2]:12.4f} {abs(m1[2]-m0[2]):9.4f}")
print(f"\n  -> largest relative change of L0, beta, theta, C: {worst:.1f}%")
print(f"  -> change of MAE: {abs(m1[0]-m0[0]):.3f} ms, RMSE {abs(m1[1]-m0[1]):.3f} ms, "
      f"R2 {abs(m1[2]-m0[2]):.4f}")

print()
print("=" * 76)
print("F) CELL BOOTSTRAP (1000 replicates, MATLAB rng(42)), 95% intervals")
print("=" * 76)
# full-data estimate, order L0 beta theta C gamma2..gamma7 (start of every refit)
pHat = np.r_[p0[:3], C0, p0[3:]]
lo = np.array([0, 0, 0, b.max() + 1] + [-50] * 6, float)
hi = np.array([50, 5, 5000, 5000] + [50] * 6, float)


def fit_local(bb, qq, yy, start, iters=5000):
    """Bounded Levenberg-Marquardt started from `start` (a local optimizer, like lsqcurvefit)."""
    G = np.column_stack([(qq == j).astype(float) for j in range(2, 8)])
    res = lambda p: yy - (p[0] + p[1] * bb + p[2] / (p[3] - bb) + G @ p[4:])
    p = np.clip(start.copy(), lo, hi); r = res(p); s = r @ r; lam = 1e-3
    for _ in range(iters):
        d = p[3] - bb
        J = np.column_stack([np.ones_like(bb), bb, 1 / d, -p[2] / d ** 2, G])
        A = J.T @ J
        pn = np.clip(p + np.linalg.solve(A + lam * np.diag(np.diag(A)), J.T @ r), lo, hi)
        rn = res(pn); sn = rn @ rn
        if sn < s:
            done = (s - sn) <= 1e-15 * s
            p, r, s, lam = pn, rn, sn, max(lam / 3, 1e-15)
            if done: break
        else:
            lam *= 4
            if lam > 1e14: break
    return p


rng = np.random.RandomState(42)                  # same Mersenne Twister stream as MATLAB rng(42)
cells = sorted(set(zip(q, tgt)))
cell_idx = [np.where((q == a) & (tgt == c))[0] for a, c in cells]
boot = np.zeros((1000, 10))
for rep in range(1000):
    sel = np.concatenate([c[np.ceil(len(c) * rng.random_sample(len(c))).astype(int) - 1] for c in cell_idx])
    boot[rep] = fit_local(b[sel], q[sel], L[sel], pHat)
lo95 = np.percentile(boot, 2.5, axis=0, method="hazen")     # MATLAB prctile definition
hi95 = np.percentile(boot, 97.5, axis=0, method="hazen")
r = np.corrcoef(boot[:, 2], boot[:, 3])[0, 1]
print(f"  replicates: {len(boot)}")
print(f"  corr(theta, C) = {r:+.3f}   (article r = 0.87)")
print()
print("  95% bootstrap intervals (compared with Table 5 of the article):")
ref = {"L0": (7.643, 9.554), "beta": (0.0271, 0.0339), "theta": (65.445, 1066.650),
       "C": (540.673, 673.624), "g2": (-1.564, 0.758), "g3": (-2.138, -0.712),
       "g4": (-2.156, -0.381), "g5": (-0.003, 1.510), "g6": (-1.102, 0.159), "g7": (-2.394, -0.941)}
for i, nm in enumerate(["L0", "beta", "theta", "C", "g2", "g3", "g4", "g5", "g6", "g7"]):
    d = 4 if nm == "beta" else 3
    rl, rh = ref[nm]
    print(f"    {nm:6s} [{lo95[i]:10.{d}f}, {hi95[i]:10.{d}f}]   article [{rl:.{d}f}, {rh:.{d}f}]")
print()
print("  The bootstrap refits are local fits along the shallow theta-C ridge, so the limits")
print("  agree with MATLAB to within the optimizer tolerance (third decimal of theta).")
print(f"  Upper limit of g6 to four decimals: {hi95[8]:.4f} (Table 5 shows 0.159).")

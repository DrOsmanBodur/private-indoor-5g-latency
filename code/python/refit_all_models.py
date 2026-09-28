"""Refits all five compared models without MATLAB and compares them with the article.

1) Full-data fit of the proposed model (Table 5, Section 5.1).
2) Trial-based train (trials 1-3) / test (trials 4-5) comparison of the five
   models (Table 6).
3) Training-set coefficients of the five models.

The nonlinearity is only in C: for a given C every model is linear in its
other parameters. The fit therefore scans C and solves a (bounded) linear
least-squares problem for each C. The bounds are those of Section 4.3.2:
0 <= L0 <= 50, 0 <= beta <= 5, 0 <= theta <= 5000, max(b) + 1 <= C <= 5000.
"""
import csv
import itertools
import os

import numpy as np

# Paths are resolved relative to this script so the package is self-contained.
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(_HERE, "..", "..", "data"))

# ---- data ----
lat = {}
for r in csv.DictReader(open(os.path.join(DATA, "latency_210.csv"), encoding="utf-8-sig")):
    lat[(int(r["qIdx"]), int(float(r["B_Mbps"])), int(r["trial"]))] = float(r["latency_ms"])

bw = {}
for r in csv.DictReader(open(os.path.join(DATA, "latency_trials_raw_210.csv"),
                             encoding="utf-8")):
    bw[(int(r["qIdx"]), int(r["target_Mbps"]), int(r["trial"]))] = float(r["achieved_Mbps"])

keys = sorted(lat)
assert len(keys) == 210 and all(k in bw for k in keys)
q = np.array([k[0] for k in keys])
tri = np.array([k[2] for k in keys])
L = np.array([lat[k] for k in keys])
b = np.array([bw[k] for k in keys])

train = np.isin(tri, [1, 2, 3])
test = np.isin(tri, [4, 5])

C_MAX = 5000.0


def lstsq(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def lstsq_bounded(X, y, lo, hi):
    """Least squares with box constraints on the first len(lo) coefficients.

    Exact for small problems: if the unconstrained solution violates a bound,
    every combination of free/lower/upper for the bounded coefficients is tried.
    """
    p = lstsq(X, y)
    k = len(lo)
    if np.all(p[:k] >= lo) and np.all(p[:k] <= hi):
        return p
    best = None
    for state in itertools.product((0, 1, 2), repeat=k):      # 0 free, 1 lower, 2 upper
        fixed = np.array([np.nan if s == 0 else (lo[i] if s == 1 else hi[i]) for i, s in enumerate(state)])
        free = np.r_[np.isnan(fixed), np.ones(X.shape[1] - k, bool)]
        pp = np.zeros(X.shape[1])
        pp[:k][~np.isnan(fixed)] = fixed[~np.isnan(fixed)]
        rhs = y - X[:, ~free] @ pp[~free]
        if free.any():
            pp[free] = lstsq(X[:, free], rhs)
        if np.all(pp[:k] >= lo - 1e-12) and np.all(pp[:k] <= hi + 1e-12):
            s = float(np.sum((y - X @ pp) ** 2))
            if best is None or s < best[0]:
                best = (s, pp)
    return best[1]


def metrics(y, yh):
    e = y - yh
    return (np.mean(np.abs(e)), np.sqrt(np.mean(e ** 2)),
            1 - np.sum(e ** 2) / np.sum((y - y.mean()) ** 2))


# ---------- design matrices (linear once C is given) ----------
def X_full(bb, qq, C):                      # L0 + beta*b + theta/(C-b) + gamma_q
    cols = [np.ones_like(bb), bb, 1.0 / (C - bb)]
    for j in range(2, 8):
        cols.append((qq == j).astype(float))
    return np.column_stack(cols)


def X_noloc(bb, qq, C):                     # L0 + beta*b + theta/(C-b)
    return np.column_stack([np.ones_like(bb), bb, 1.0 / (C - bb)])


def X_mm1(bb, qq, C):                       # L0 + theta/(C-b)
    return np.column_stack([np.ones_like(bb), 1.0 / (C - bb)])


def X_poly(bb, qq, C=None):                 # poly22
    return np.column_stack([np.ones_like(bb), bb, qq, bb ** 2, bb * qq, qq ** 2])


def X_persat(bb, qq, C=None, Bth=200.0):    # L0(q)+alpha(q)ln(1+b)+beta(q)max(0,b-Bth)
    cols = []
    f1, f2 = np.log1p(bb), np.maximum(0.0, bb - Bth)
    for j in range(1, 8):
        m = (qq == j).astype(float)
        cols += [m, m * f1, m * f2]
    return np.column_stack(cols)


# bounds of the leading coefficients of each C-model (Section 4.3.2)
BOUNDS = {X_full: ([0, 0, 0], [50, 5, 5000]),
          X_noloc: ([0, 0, 0], [50, 5, 5000]),
          X_mm1: ([0, 0], [50, 5000])}


def fit_profiled(Xfun, bb, qq, yy):
    """1-D scan over C with a bounded linear solve for every C, then local refinement."""
    lo_b, hi_b = (np.array(v, float) for v in BOUNDS[Xfun])
    lo = bb.max() + 1.0

    def sse(C):
        X = Xfun(bb, qq, C)
        p = lstsq_bounded(X, yy, lo_b, hi_b)
        r = yy - X @ p
        return float(r @ r), p

    best = None
    for C in np.concatenate([np.linspace(lo, lo + 50, 400), np.linspace(lo + 50, C_MAX, 3000)]):
        s, p = sse(C)
        if best is None or s < best[0]:
            best = (s, C, p)
    s, C, p = best
    step = 1.0
    for _ in range(60):
        improved = False
        for cand in (C - step, C + step):
            if cand < lo or cand > C_MAX:
                continue
            ss, pp = sse(cand)
            if ss < s:
                s, C, p, improved = ss, cand, pp, True
        if not improved:
            step /= 2
            if step < 1e-6:
                break
    return C, p


MODELS = {
    "(a) proposed full":       (X_full,   True),
    "(b) no-location":         (X_noloc,  True),
    "(c) minimal M/M/1-type":  (X_mm1,    True),
    "(d) polynomial poly22":   (X_poly,   False),
    "(e) per-location sat.":   (X_persat, False),
}

print("=" * 78)
print("1) FULL DATASET, MODEL (a)  --  compared with Table 5 of the article")
print("=" * 78)
C, p = fit_profiled(X_full, b, q, L)
X = X_full(b, q, C)
mae, rmse, r2 = metrics(L, X @ p)
print(f"  L0    = {p[0]:8.3f}  ms        (article 8.639)")
print(f"  beta  = {p[1]:8.4f}  ms/Mbps   (article 0.0305)")
print(f"  theta = {p[2]:8.3f}  ms*Mbps   (article 486.885)")
print(f"  C     = {C:8.3f}  Mbps      (article 611.395)")
gp = [-0.470, -1.449, -1.306, 0.782, -0.476, -1.681]
for i, j in enumerate(range(2, 8)):
    print(f"  gamma{j} = {p[3 + i]:8.3f}  ms        (article {gp[i]:+.3f})")
print(f"  MAE={mae:.3f}  RMSE={rmse:.3f}  R2={r2:.4f}   (article 2.52 / 3.13 / 0.842)")

print()
print("=" * 78)
print("2) TRIAL-BASED TRAIN (1-3) / TEST (4-5)  --  compared with Table 6 of the article")
print("=" * 78)
ref = {"(a) proposed full": (2.654, 3.269, 0.835),
       "(b) no-location": (2.655, 3.329, 0.829),
       "(c) minimal M/M/1-type": (2.892, 3.633, 0.797),
       "(d) polynomial poly22": (2.724, 3.382, 0.824),
       "(e) per-location sat.": (1.787, 2.322, 0.917)}
coef = {}
for name, (Xf, hasC) in MODELS.items():
    if hasC:
        Ci, p = fit_profiled(Xf, b[train], q[train], L[train])
    else:
        Ci = None
        p = lstsq(Xf(b[train], q[train]), L[train])
    yh = (Xf(b[test], q[test], Ci) if hasC else Xf(b[test], q[test])) @ p
    m = metrics(L[test], yh)
    coef[name] = (Ci, p)
    r = ref[name]
    print(f"  {name:24s} test MAE={m[0]:6.3f} RMSE={m[1]:6.3f} R2={m[2]:6.3f}"
          f"   | article {r[0]:.3f}/{r[1]:.3f}/{r[2]:.3f}")

print()
print("=" * 78)
print("3) TRAINING-SET COEFFICIENTS")
print("=" * 78)
Ci, p = coef["(a) proposed full"]
print("(a) proposed full model")
print(f"    L0={p[0]:.3f}  beta={p[1]:.4f}  theta={p[2]:.2f}  C={Ci:.2f}")
print("    " + "  ".join(f"g{j}={p[3+i]:+.3f}" for i, j in enumerate(range(2, 8))))
Ci, p = coef["(b) no-location"]
print(f"(b) no-location ablation\n    L0={p[0]:.3f}  beta={p[1]:.4f}  theta={p[2]:.2f}  C={Ci:.2f}")
Ci, p = coef["(c) minimal M/M/1-type"]
print(f"(c) minimal M/M/1-type (theta at its upper bound of 5000)\n    L0={p[0]:.3f}  theta={p[1]:.2f}  C={Ci:.2f}")
_, p = coef["(d) polynomial poly22"]
print("(d) polynomial poly22")
print("    " + "  ".join(f"a{i}={v:+.5g}" for i, v in enumerate(p)))
_, p = coef["(e) per-location sat."]
print("(e) per-location saturation (Bth=200 Mbps)")
print("    q :      L0(q)    alpha(q)   beta(q)")
for j in range(7):
    print(f"    Q{j+1}: {p[3*j]:9.3f} {p[3*j+1]:10.3f} {p[3*j+2]:9.4f}")

"""Chowla cosine problem tools. Score = -min f_A / sqrt(|A|) (lower is better).

  python chowla.py baseline
  python chowla.py search --n 200 --L 4000 --iters 50000 [--seed 0]
"""
import argparse, json, math, os, random, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

def fmin(A, oversample=16):
    """Min of f_A over [0, pi]: FFT grid, then Newton refinement at the lowest grid points."""
    A = np.asarray(sorted(A), dtype=np.int64)
    M = 1 << int(math.ceil(math.log2(oversample * (A.max() + 1))))
    c = np.zeros(M); c[A] = 1.0
    f = np.fft.rfft(c).real            # f[k] = f_A(2*pi*k/M)
    idx = np.argsort(f)[:8]
    best = f[idx[0]]
    for k in idx:
        x = 2 * math.pi * k / M
        for _ in range(20):
            s = np.sin(A * x); co = np.cos(A * x)
            d1 = -(A * s).sum(); d2 = -(A * A * co).sum()
            if d2 <= 0: break
            x -= d1 / d2
        best = min(best, np.cos(A * x).sum())
    return float(best)

def score(A):
    return -fmin(A) / math.sqrt(len(A))

def is_prime(p):
    return p > 1 and all(p % q for q in range(2, int(p ** .5) + 1))

def baseline():
    rows = []
    for p in [101, 401, 1009, 4001]:
        if not is_prime(p): continue
        qr = sorted({x * x % p for x in range(1, p)})
        rows.append((f"QR mod {p}", len(qr), score(qr)))
    for n in [100, 400, 1000]:
        rows.append((f"squares 1..{n}^2", n, score([i * i for i in range(1, n + 1)])))
        rows.append((f"interval 1..{n}", n, score(range(1, n + 1))))
        rnd = random.Random(0).sample(range(1, 4 * n), n)
        rows.append((f"random n={n} in [1,{4*n})", n, score(rnd)))
    for name, n, s in rows:
        print(f"{name:28s} |A|={n:5d}  score={s:.4f}")

def search(n, L, iters, seed, init=None):
    rng = random.Random(seed)
    A = set(init) if init else set(rng.sample(range(1, L + 1), n))
    n = len(A)
    cur = score(A); best, bestA = cur, set(A)
    T0, t0 = (0.01 if init else 0.05), time.time()
    for it in range(iters):
        T = T0 * (1 - it / iters) + 1e-4
        B = set(A)
        out = rng.choice(tuple(B)); B.remove(out)
        while True:
            inn = rng.randint(1, L)
            if inn not in B: B.add(inn); break
        s = score(B)
        if s < cur or rng.random() < math.exp((cur - s) / T):
            A, cur = B, s
            if s < best:
                best, bestA = s, set(B)
        if it % 2000 == 0:
            print(f"it={it} cur={cur:.4f} best={best:.4f} ({time.time()-t0:.0f}s)", flush=True)
    os.makedirs(os.path.join(HERE, "runs"), exist_ok=True)
    tag = f"{'pip' if init else ''}n{n}_L{L}_s{seed}"
    with open(os.path.join(HERE, "runs", tag + ".json"), "w") as fh:
        json.dump({"n": n, "L": L, "score": best, "A": sorted(bestA)}, fh)
    with open(os.path.join(HERE, "runs", "best.txt"), "a") as fh:
        fh.write(f"{time.strftime('%F %T')} {tag} score={best:.4f}\n")
    print("BEST", tag, best)

def pipage_round(K, beta=8.0, seed=0, oversample=8):
    """Start from Fejer weights w_k = 1 - k/K on 1..K (min f >= -1/2), round to 0/1 by pipage
    moves along e_i - e_j that keep sum fixed and push up a soft-min of f over a grid."""
    rng = np.random.default_rng(seed)
    k = np.arange(1, K + 1)
    w = 1.0 - k / K
    M = 1 << int(math.ceil(math.log2(oversample * (K + 1))))
    xs = 2 * np.pi * np.arange(M // 2 + 1) / M
    C = np.cos(np.outer(xs, k))                      # grid x element
    frac = lambda: np.where((w > 1e-9) & (w < 1 - 1e-9))[0]
    F = frac()
    while len(F) >= 2:
        f = C @ w
        n = w.sum()
        z = -beta * f / math.sqrt(n); z -= z.max(); p = np.exp(z); p /= p.sum()
        g = p @ C                                     # gradient of soft-min wrt w
        i, j = rng.choice(F, 2, replace=False)
        d = 1.0 if g[i] > g[j] else -1.0              # move along d*(e_i - e_j)
        step = min(1 - w[i], w[j]) if d > 0 else min(w[i], 1 - w[j])
        w[i] += d * step; w[j] -= d * step
        w[np.abs(w) < 1e-9] = 0; w[np.abs(w - 1) < 1e-9] = 1
        F = frac()
    if len(F) == 1: w[F[0]] = 1.0 if rng.random() < w[F[0]] else 0.0
    return [int(a) for a in k[w > 0.5]]

def rounding(Ks, seeds):
    for K in Ks:
        best = min((score(pipage_round(K, seed=s)), s) for s in range(seeds))
        A = pipage_round(K, seed=best[1])
        print(f"K={K:5d} |A|={len(A):5d} score={best[0]:.4f}", flush=True)
        with open(os.path.join(HERE, "runs", "best.txt"), "a") as fh:
            fh.write(f"{time.strftime('%F %T')} pipage K{K} n{len(A)} score={best[0]:.4f}\n")

def lm_round(K, tc=0.7, gamma=0.03, eta=1.0, seed=0, oversample=8, verbose=False):
    """Lovett-Meka style rounding with a drift. Start at Fejer weights; take Gaussian + soft-min
    gradient steps projected orthogonal to (a) the sum constraint and (b) grid constraints that are
    near the floor -tc*sqrt(n). Coordinates freeze at 0/1. Leftovers finished by pipage."""
    rng = np.random.default_rng(seed)
    k = np.arange(1, K + 1); w = 1.0 - k / K
    M = 1 << int(math.ceil(math.log2(oversample * (K + 1))))
    xs = 2 * np.pi * np.arange(M // 2 + 1) / M
    C = np.cos(np.outer(xs, k))
    n0 = w.sum(); floor = -tc * math.sqrt(n0)
    for it in range(200000):
        F = np.where((w > 1e-6) & (w < 1 - 1e-6))[0]
        if verbose and it % 200 == 0: print(it, len(F), w[F].min(), w[F].max(), flush=True)
        if len(F) < 4: break
        f = C @ w
        T = np.where(f < floor + 0.15 * math.sqrt(n0))[0]
        T = T[np.argsort(f[T])][: max(0, len(F) // 2 - 1)]
        B = np.vstack([np.ones(len(F)), C[np.ix_(T, F)]]) if len(T) else np.ones((1, len(F)))
        Q, _ = np.linalg.qr(B.T)
        z = -8 * f / math.sqrt(n0); z -= z.max(); p = np.exp(z); p /= p.sum()
        d = rng.standard_normal(len(F)) + eta * math.sqrt(len(F)) * (p @ C[:, F]) / (np.linalg.norm(p @ C[:, F]) + 1e-12)
        d -= Q @ (Q.T @ d)
        nd = np.linalg.norm(d)
        if nd < 1e-9: break
        d *= gamma / nd * math.sqrt(len(F))
        lim = np.where(d > 0, (1 - w[F]) / np.maximum(d, 1e-12), w[F] / np.maximum(-d, 1e-12))
        s = min(1.0, lim.min())
        if verbose and it % 200 == 0: print("  s", s, np.abs(d).max())
        w[F] += s * d
        w = np.clip(w, 0, 1); w[w < 1e-6] = 0; w[w > 1 - 1e-6] = 1
    # finish with pipage on the leftovers
    F = np.where((w > 0) & (w < 1))[0]
    while len(F) >= 2:
        f = C @ w; z = -8 * f / math.sqrt(n0); z -= z.max(); p = np.exp(z); p /= p.sum(); g = p @ C
        i, j = rng.choice(F, 2, replace=False); dd = 1.0 if g[i] > g[j] else -1.0
        st = min(1 - w[i], w[j]) if dd > 0 else min(w[i], 1 - w[j])
        w[i] += dd * st; w[j] -= dd * st
        w[np.abs(w) < 1e-9] = 0; w[np.abs(w - 1) < 1e-9] = 1
        F = np.where((w > 0) & (w < 1))[0]
    if len(F) == 1: w[F[0]] = 1.0
    return [int(a) for a in k[w > 0.5]]

def lm_sweep(Ks, tcs, seeds):
    for K in Ks:
        for tc in tcs:
            best = min((score(lm_round(K, tc=tc, seed=s)), s) for s in range(seeds))
            line = f"lm K{K} tc{tc} score={best[0]:.4f} seed={best[1]}"
            print(line, flush=True)
            with open(os.path.join(HERE, "runs", "best.txt"), "a") as fh:
                fh.write(f"{time.strftime('%F %T')} {line}\n")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["baseline", "search", "round", "polish", "lm"])
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--L", type=int, default=2000)
    ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "baseline": baseline()
    elif a.cmd == "round": rounding([200, 400, 1000, 2000], 4)
    elif a.cmd == "lm":
        lm_sweep([a.L], [0.4, 0.6, 0.8], 3)
    elif a.cmd == "polish":
        A = pipage_round(a.L, seed=a.seed); search(len(A), a.L, a.iters, a.seed, init=A)
    else: search(a.n, a.L, a.iters, a.seed)

"""Continuous melt-down: optimize real weights w_a in [0,1] (a = 1..L, sum w = n) to maximize a soft-min
of f(x) = sum w_a cos(a x), while a 0/1 penalty mu*sum w(1-w) ramps up so w crystallizes into a set.
f and gradients via FFT, O(G log G) per step.   python melt.py n L SECONDS SEED [mumax] [tau_end]"""
import json, math, os, sys, time
import numpy as np
import chowla as c

def run(n, L, seconds, seed, mumax=3.0, tau_end=0.3, init=None):
    rng = np.random.default_rng(seed); t0 = time.time()
    if init: L = max(L, max(json.load(open(init))['A']))
    G = 1 << int(math.ceil(math.log2(8 * L))); a = np.arange(1, L + 1)
    N = 2 * n
    if init:
        A0 = json.load(open(init))['A']; L = max(L, max(A0)); a = np.arange(1, L + 1); G = 1 << int(math.ceil(math.log2(8 * L)))
        w = 0.5 + 0.3 * (np.isin(a, A0) - 0.5) + 0.05 * rng.standard_normal(L)
    else:
        w = np.clip(1 - a / N, 0, 1) + 0.05 * rng.standard_normal(L)
    def proj(w):
        lo, hi = -1.0, 1.0
        for _ in range(50):
            s = (lo + hi) / 2
            if np.clip(w + s, 0, 1).sum() > n: hi = s
            else: lo = s
        return np.clip(w + lo, 0, 1)
    w = proj(w); buf = np.zeros(2 * G)
    def fvals(w):
        buf[:] = 0; buf[1:L + 1] = w; return np.fft.rfft(buf).real[:G + 1]
    it = 0; best = 9e9; bestA = None
    while True:
        frac = (time.time() - t0) / seconds
        if frac >= 1: break
        it += 1
        tau = 3.0 * (tau_end / 3.0) ** frac; mu = mumax * frac ** 2
        f = fvals(w); fm = f.min()
        p = np.exp(-(f - fm) / tau); p /= p.sum()
        q = np.zeros(2 * G); q[:G + 1] = p; q[G + 1:] = p[1:G][::-1]
        g = np.fft.rfft(q).real[1:L + 1] / 2          # d softmin / d w_a  ~ sum_j p_j cos(a x_j)
        g = g / (np.abs(g).max() + 1e-12) + mu * (2 * w - 1)
        w = proj(w + 0.02 * g)
        if it % 200 == 0:                              # snapshot: round to top-n and score
            A = (np.argsort(-w)[:n] + 1).tolist(); fa = fvals(np.isin(a, A).astype(float)).min()
            sc = -fa / math.sqrt(n)
            if sc < best: best, bestA = sc, A
    A = sorted(bestA); sc = c.score(A); dense = -c.fmin(A, oversample=64) / math.sqrt(n)
    line = f"MELT init={os.path.basename(init) if init else 'fejer'} n={n} L={L} mumax={mumax} tau_end={tau_end} score={sc:.4f} dense={dense:.4f} iters={it}"
    print(line, flush=True)
    json.dump({"A": A, "dense": dense}, open(os.path.join(c.HERE, "runs", f"melt_n{n}_L{L}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh: fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    v = sys.argv[1:]
    run(int(v[0]), int(v[1]), float(v[2]), int(v[3]), float(v[4]), float(v[5]), v[6] if len(v) > 6 and v[6] != '-' else None)

"""COMBO: fixed optimized Sidon core + LP correction layer + one-sided pipage rounding + swap polish.
  python correct.py CORE.json MFRAC SECONDS SEED
Score = -min f / sqrt(|A|) over the final set A = core U E (E disjoint from core)."""
import json, math, os, sys, time
import numpy as np
from scipy.optimize import linprog
import chowla as c

def main(core_file, mfrac, seconds, seed):
    rng = np.random.default_rng(seed); t0 = time.time()
    core = sorted(json.load(open(core_file))["A"]); n0 = len(core); L = max(core)
    m = max(1, int(round(mfrac * n0)))
    G = 1 << int(math.ceil(math.log2(4 * L))); x = np.pi * np.arange(G + 1) / G
    ind = np.zeros(2 * G); ind[core] = 1
    fcore = np.fft.rfft(ind).real[:G + 1]                  # f_core on x_j = pi*j/G
    cand = np.setdiff1d(np.arange(1, L + 1), core)         # correction frequencies
    # cutting-plane LP: max t s.t. fcore(x_j) + sum_e w_e cos(e x_j) >= t on active points, sum w = m
    act = np.argsort(fcore)[:1500]; w = None
    rnd = -1
    while time.time() - t0 < 0.7 * seconds:
        rnd += 1
        C = np.cos(np.outer(x[act], cand))
        nv = len(cand)
        A_ub = np.hstack([-C, np.ones((len(act), 1))]); b_ub = fcore[act]
        res = linprog(np.r_[np.zeros(nv), -1.0], A_ub=A_ub, b_ub=b_ub,
                      A_eq=np.r_[np.ones(nv), 0.0][None, :], b_eq=[m],
                      bounds=[(0, 1)] * nv + [(None, None)], method="highs")
        w = res.x[:nv]
        wf = np.zeros(2 * G); wf[cand] = w
        ftot = fcore + np.fft.rfft(wf).real[:G + 1]
        new = np.setdiff1d(np.argsort(ftot)[:400], act)
        print(f"LP round {rnd}: t={res.x[-1]:.2f} grid-min={ftot.min():.2f} frac-score={-ftot.min()/math.sqrt(n0+m):.4f} ({time.time()-t0:.0f}s)", flush=True)
        if len(new) == 0: break
        act = np.r_[act, new]
    # one-sided pipage rounding on fractional coords, O(G) per step
    F = np.flatnonzero((w > 1e-6) & (w < 1 - 1e-6))
    wf = np.zeros(2 * G); wf[cand] = w; f = fcore + np.fft.rfft(wf).real[:G + 1]
    soft = lambda v: -np.log(np.exp(-(v - v.min()) / 2.0).sum()) * 2.0 + v.min()
    w = w.copy(); fr = list(F); rng.shuffle(fr)
    while len(fr) >= 2:
        i, j = fr[-1], fr[-2]
        ci, cj = np.cos(cand[i] * x), np.cos(cand[j] * x)
        up = min(1 - w[i], w[j]); dn = min(w[i], 1 - w[j])
        fa = f + up * (ci - cj); fb = f - dn * (ci - cj)
        if soft(fa) >= soft(fb): w[i] += up; w[j] -= up; f = fa
        else: w[i] -= dn; w[j] += dn; f = fb
        fr = [k for k in fr if 1e-9 < w[k] < 1 - 1e-9]
    E = cand[w > 0.5]
    A = sorted(set(core) | set(int(e) for e in E))
    # polish: swap one correction element for an unused frequency, keep if min improves
    Aset = set(A); Ecur = [int(e) for e in E]
    ind = np.zeros(2 * G); ind[A] = 1; f = np.fft.rfft(ind).real[:G + 1]; cur = f.min()
    while time.time() - t0 < seconds and Ecur:
        a = Ecur[rng.integers(len(Ecur))]; b = int(cand[rng.integers(len(cand))])
        if b in Aset: continue
        f2 = f - np.cos(a * x) + np.cos(b * x); v = f2.min()
        if v >= cur:
            f, cur = f2, v; Aset.remove(a); Aset.add(b); Ecur[Ecur.index(a)] = b
    A = sorted(Aset); n = len(A)
    sc = c.score(A); dense = -c.fmin(A, oversample=64) / math.sqrt(n)
    tag = os.path.basename(core_file).replace(".json", "")
    line = f"CORR {tag} mfrac={mfrac} n0={n0} n={n} score={sc:.4f} dense={dense:.4f}"
    print(line, flush=True)
    json.dump({"core": core_file, "mfrac": mfrac, "A": A, "score": sc, "dense": dense},
              open(os.path.join(c.HERE, "runs", f"corr_{tag}_m{mfrac}_s{seed}.json"), "w"))
    with open(os.path.join(c.HERE, "runs", "best.txt"), "a") as fh:
        fh.write(time.strftime("%F %T ") + line + "\n")

if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))

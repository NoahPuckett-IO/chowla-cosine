# Notes from Mac Claude (2026-10-17)
Read NOTES_FROM_PC.md. Agreed: BLAS-level chaos. Here Accelerate gave 0.906 for the K=300 test and OpenBLAS gives 0.882; same code.

## Mac setup
- ~/math/venv: uv-managed CPython 3.12, numpy 2.5.4 **OpenBLAS wheel**. The default macOS wheel uses Accelerate and corrupted the heap here (pthread_mutex errors, SIGTRAP in malloc).
- Jobs run at nice 20 with OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1. Scripts are unchanged.
- Only 8 GiB RAM, so the big K values run one at a time and K=10000 is skipped (its C matrix alone is ~4.9 GiB).

## Results so far (Mac, sweep_lm.py)
| K | n | score | (eta, tc, seed) |
|---|---|---|---|
| 1200 | 600 | **0.7783** | (0.5, 0.4, 2) |
| 1000 | 500 | 0.7969 | (0.5, 0.5, 0) |
| 1400 | 700 | 0.8023 | (0.2, 0.4, 2) |

## Running / queued (as of 13:30)
- K=5000 → 6000 → 8000 (one at a time; 5000 started 11:58)
- K=2000 → 2500 → 3000 → 4000 (one at a time; 2000 started 12:07)
- K=1700 (started 13:23)

If you're running the same K values on the PC, note it here so we don't duplicate work. I'll keep this file updated.

## Reply to PC ~13:40 update (14:10)
- Got it. Keeping 3000/4000 here: only 3 of 6 Mac cores are busy, so cores aren't tight, and the duplicates give extra samples from a different BLAS.
- Not duplicated by the PC: K=1700 and the 5000 → 6000 → 8000 chain.
- Mac awake: `caffeinate -is -w` is tied to both long chains, so it lasts until they finish.
- Status 14:06: K=1700 at 42m, K=2000 at 1h58m, K=5000 at 2h08m. Nothing new finished since 13:23.

## 14:12 — all 6 Mac cores in use (Noah's request). PC Claude is primary; redirect me here if you want different K values.
Added chains, chosen to avoid your K values and kept at K ≤ 3500 to fit in RAM:
- 1100 → 1900 → 2800
- 1300 → 2200 → 3500
- 1500 → 900 → 2400

## Results log (Mac)
- 14:12 K2000 n1000 score=0.7839 params=(0.3, 0.3, 1)
- 15:13 K1100 n550 score=0.7672 params=(0.3, 0.3, 0)

## 15:20: all Mac Chowla sweeps STOPPED by Noah; they will not be restarted
~0.78 is beaten by plain Sidon differences (~0.707) and by Bu et al. (arXiv 2610.17160, constant 1/2). Summary, table and construction: BU_ET_AL.md.
Side note: around 15:13 every background job here was killed at once, logs clean, cause unknown.

## BUOPT results (Mac, 40 min each, seed 2; finished 16:04)
- 2026-10-17 16:03:56 BUOPT k40 p53 n1560 grid=0.6518 score=0.6519 dense=0.6519 iters=2661661
- 2026-10-17 16:03:57 BUOPT k40 p71 n1560 grid=0.6618 score=0.6619 dense=0.6619 iters=1279400
- 2026-10-17 16:03:57 BUOPT k60 p73 n3540 grid=0.6269 score=0.6282 dense=0.6282 iters=1275567
- 2026-10-17 16:03:57 BUOPT k60 p79 n3540 grid=0.6250 score=0.6271 dense=0.6271 iters=1273599
- 2026-10-17 16:03:59 BUOPT k80 p89 n6320 grid=0.6143 score=0.6167 dense=0.6167 iters=602124
- 2026-10-17 16:04:00 BUOPT k80 p101 n6320 grid=0.6142 score=0.6171 dense=0.6171 iters=597985
- Independent check on all 6 saved sets: S is Sidon, A has k(k−1) distinct positive elements, and min f_A recomputed by 64× FFT plus direct refinement matches the saved dense score.
- vs Bu Table 1: n1560 0.6519 (vs 0.725), n3540 0.6271 (vs 0.708), n6320 0.6167 (vs 0.674). All beaten. Best p per k: 53 / 79 / 89.
- All 6 cores are free now.

## BUOPT round 2 (Mac, 60 min each; finished 17:07)
- 2026-10-17 17:06:10 BUOPT k100 p127 n9900 grid=0.6016 score=0.6050 dense=0.6050 iters=288128
- 2026-10-17 17:06:38 BUOPT k100 p113 n9900 grid=0.6043 score=0.6058 dense=0.6058 iters=681233
- 2026-10-17 17:07:06 BUOPT k40 p47 n1560 grid=0.6433 score=0.6450 dense=0.6450 iters=3901433
- 2026-10-17 17:07:07 BUOPT k60 p67 n3540 grid=0.6242 score=0.6242 dense=0.6242 iters=1770414
- 2026-10-17 17:07:13 BUOPT k80 p89 n6320 grid=0.6138 score=0.6140 dense=0.6140 iters=711967
- 2026-10-17 17:07:13 BUOPT k80 p97 n6320 grid=0.6162 score=0.6169 dense=0.6169 iters=697144
- Independent check on all 6 round-2 sets: PASS (Sidon S, k(k−1) distinct positive elements, recomputed min matches saved dense score).
- All 6 Mac cores are free now.

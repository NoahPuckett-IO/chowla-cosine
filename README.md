# Chowla's cosine problem - workspace index

Problem: Epoch FrontierMath Open Problem. Score = -min f_A / sqrt(|A|), lower is better. See PROBLEM.md.

## Results (final table, 2026-10-17, dense-verified on PC + Mac)
| n | Ours (optimized Sidon selection in Bu identity) | Bu et al. 2026 (random selection) |
|---|---|---|
| 1,560 | 0.645 | 0.725 |
| 2,450 | 0.640 | - |
| 3,540 | 0.624 | 0.708 |
| 4,830 | 0.621 | - |
| 6,320 | 0.614 | 0.674 |
| 8,010 | 0.607 | - |
| 9,900 | 0.605 | - |
| 11,990 | 0.598 | - |
| 14,280 | 0.598 | 0.655 |
| 25,440 | 0.590 | 0.635 |
Classical Sidon differences: ~0.71. Bu et al. prove the limit 1/2; our method approaches it from above.

## Files
- FINDINGS.md - method + all result tables. LITERATURE.md - what is known (Bu et al. arXiv 2610.17160).
- IDEAS.md - ranked list of 7 ideas toward much lower scores / disproof.
- PROMPT_FOR_PC.md - morning prompt for PC Claude.
- NEXT_STRATEGY.md - proposed "Sidon core + correction layer" attack on the 1/2 floor.
- PLAN.md - day log. PROBLEM.md - statement.
- bu_opt.py - current best method. chowla.py - scorer + earlier methods (pipage, LM rounding).
- sweep_lm.py - old LM sweeps (superseded). PROMPT_FOR_MAC.md - task file mirrored to mac:~/math.
- runs/best.txt, runs/buopt_*.json (PC sets); runs/mac/ (Mac results, sets, notes).
Mac: ssh mac, code in ~/math/Chowla, Python at ~/math/venv (OpenBLAS numpy).

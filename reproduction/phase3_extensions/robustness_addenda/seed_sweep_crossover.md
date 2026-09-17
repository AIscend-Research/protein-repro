# Seed-stability of the noise-sweep crossover

Each cell is the mean over 3 seeds (37, 38, 39) x 30 structures x 8 samples at T=0.1 (240 sequences per seed); noise=0.5/seed=37 reuses the original single-seed noise_sweep run rather than re-computing it.

| noise (A) | full-backbone mean (range across seeds) | CA-only mean (range across seeds) | FB - CA (pp) | which model wins |
|---|---|---|---|---|
| 0.4 | 32.3% (32.2-32.5%) | 30.4% (30.2-30.6%) | +1.9 | full-backbone |
| 0.45 | 28.3% (28.3-28.4%) | 27.9% (27.8-28.2%) | +0.4 | full-backbone |
| 0.5 | 24.3% (23.9-24.7%) | 25.3% (25.1-25.4%) | -1.0 | CA-only |
| 0.6 | 17.8% (17.7-18.0%) | 20.6% (20.3-21.0%) | -2.8 | CA-only |

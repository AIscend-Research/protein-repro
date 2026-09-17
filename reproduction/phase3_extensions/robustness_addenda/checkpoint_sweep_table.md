# CA-only recovery across training-noise checkpoints

v_48_020/v_48_002/v_48_010 = CA-only checkpoints trained with Gaussian coordinate noise std. 0.20/0.02/0.10 A respectively. v_48_020 column reuses the original noise_sweep (single seed=37) run; v_48_002/v_48_010 use seed=37 only (this sweep's purpose is to see whether checkpoint choice shifts the curve, not to re-establish seed stability, which the seed_sweep already covers).

| noise (A) | v_48_020 (matched, seed=37) | v_48_002 | v_48_010 |
|---|---|---|---|
| 0.0 | 40.5% | 42.9% | 42.5% |
| 0.2 | 38.0% | 25.0% | 36.4% |
| 0.3 | 35.1% | 19.0% | 29.2% |
| 0.4 | n/a | 15.7% | 23.5% |
| 0.5 | 25.4% | 13.7% | 19.8% |

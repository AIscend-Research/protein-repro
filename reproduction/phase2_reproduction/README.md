# Phase 2: Reproduction

Reproduces the main ProteinMPNN inference pipeline on the Phase 1 dataset (30
PDBs, see [reproduction/phase1_data_prep/README.md](../phase1_data_prep/README.md)):
sequence recovery, model score / perplexity, and runtime for the full-backbone
model across a temperature sweep, compared against the CA-only model.

## Commands

Chain assignments (`chain_id.jsonl`) mark 1TUP's DNA chains E/F as fixed and
all other chains across all 30 structures as designed (`build_chain_id_jsonl.py`).

Full-backbone model:
```
python protein_mpnn_run.py \
  --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
  --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
  --out_folder reproduction/phase2_reproduction/fullbackbone \
  --model_name v_48_020 \
  --num_seq_per_target 8 \
  --sampling_temp "0.1 0.2 0.3" \
  --seed 37 \
  --batch_size 1 \
  --save_score 1
```

CA-only model (identical command plus `--ca_only`, same `--jsonl_path` —
`tied_featurize` subsets the already-parsed full-atom backbone down to CA
coordinates internally; no separate CA-only parsing step exists or is needed):
```
python protein_mpnn_run.py --ca_only \
  --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
  --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
  --out_folder reproduction/phase2_reproduction/ca_only \
  --model_name v_48_020 \
  --num_seq_per_target 8 \
  --sampling_temp "0.1 0.2 0.3" \
  --seed 37 \
  --batch_size 1 \
  --save_score 1
```

Both runs: 30 proteins x 3 temperatures x 8 samples/target = 720 designed
sequences per model, 1440 total. Model checkpoint `v_48_020` for both
(48 edges, 0.20 A training noise) so the full-backbone vs. CA-only comparison
isn't confounded by a different noise level. Git commit
`8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`.

## Metrics pipeline
- `compute_metrics.py` parses every `seqs/*.fa` header
  (`score`, `global_score`, `seq_recovery`) into `per_sample_metrics.csv`
  (1440 rows), computes perplexity as `exp(score)` (score = average per-residue
  negative log-probability in nats, confirmed against `.npz` score files in
  Phase 0/1), and parses the terminal logs' `N sequences ... generated in T
  seconds` lines into `runtime_metrics.csv` (60 rows, one per protein per model).
- `temperature_sweep_summary.csv` aggregates mean seq_recovery / score /
  perplexity per (model, temperature).
- `make_plots.py` renders `temperature_sweep_comparison.png` (recovery and
  perplexity vs. temperature, both models) and `runtime_comparison.png`.

## Results

| model | temperature | mean seq_recovery | mean score (NLL) | mean perplexity |
|-------|-------------|--------------------|--------------------|-------------------|
| fullbackbone | 0.1 | 46.4% | 0.866 | 2.397 |
| fullbackbone | 0.2 | 45.9% | 0.898 | 2.478 |
| fullbackbone | 0.3 | 45.3% | 0.954 | 2.629 |
| ca_only | 0.1 | 41.3% | 0.953 | 2.623 |
| ca_only | 0.2 | 40.5% | 0.986 | 2.714 |
| ca_only | 0.3 | 39.9% | 1.051 | 2.909 |

**Qualitative findings reproduced:**
1. Full-backbone recovery > CA-only recovery at every temperature (+5.1 pp at
   T=0.1), consistent with the paper's claim that full backbone context
   (N, CA, C, O) improves design accuracy over CA-only.
2. Recovery decreases and perplexity increases monotonically with temperature
   for both models — the expected diversity/accuracy tradeoff.
3. Absolute recovery (~40-46%) is lower than the paper's headline ~52% on its
   large-scale CATH test set — expected, since our 30-protein set is small,
   hand-curated from RCSB (not the paper's held-out CATH-clustered test split,
   which isn't separately downloadable — see `training/README.md`, the full
   16.5 GB training tarball is the only source of `test_clusters.txt`), and
   includes some short/simple domains (e.g. 1CRN, 46 residues) that skew
   differently than the paper's length distribution.

**Input completeness (clean vs. stress-test subset) — added on re-verification,
was missing from the first pass.** `completeness_comparison.py` joins
`per_sample_metrics.csv` against the Phase 1 `dataset_summary.csv` subset
labels (21 clean / 9 stress-test proteins) and aggregates recovery/score by
(model, temperature, subset). Result is **not a clean, predictable trend**:

| model | temp | clean recovery | stress-test recovery | delta |
|---|---|---|---|---|
| fullbackbone | 0.1 | 47.0% | 44.9% | -2.1 pp |
| fullbackbone | 0.2 | 46.6% | 44.4% | -2.2 pp |
| fullbackbone | 0.3 | 45.8% | 44.1% | -1.7 pp |
| ca_only | 0.1 | 40.7% | 42.6% | **+2.0 pp** |
| ca_only | 0.2 | 39.8% | 42.0% | **+2.1 pp** |
| ca_only | 0.3 | 39.7% | 40.3% | **+0.5 pp** |

Full-backbone shows the naively-expected direction (recovery slightly lower
on stress-test structures) but the gap is small (~2 pp). CA-only shows the
**opposite** direction — stress-test structures score *better*. See
`completeness_comparison.png`.

This is a genuine finding, not a bug, but it means the roadmap's expected
"model behavior changes predictably with input completeness" does **not**
hold up under this natural-subset comparison, for two identifiable reasons:
1. ProteinMPNN's own `score`/`seq_recovery` are computed only over resolved
   (mask=1) positions — unresolved/gap residues are excluded from the loss
   entirely (`_scores()`, `protein_mpnn_utils.py:39`), so having *some*
   missing residues elsewhere in a structure doesn't directly penalize the
   score on the positions that *are* evaluated.
2. The clean/stress split confounds completeness with protein identity: 21
   vs. 9 different proteins, differing in fold, size, and composition,
   likely swamp any true completeness effect at this sample size.

A proper test of this claim needs a **matched, controlled** comparison —
same backbone, varying degrees of induced incompleteness — which is exactly
what Phase 3's synthetic residue-masking/noise experiments are designed to
do. Flagging this here so it isn't silently dropped: the natural-subset
comparison alone is not sufficient evidence either way, and Phase 3 is where
this claim should actually be tested.

**Runtime** (CPU only, Apple M3, 16GB RAM — the official code's device
selection, `protein_mpnn_run.py:68`, checks only `torch.cuda.is_available()`
and never MPS, so this Mac's GPU went unused; a deliberate low-resource
condition either way, see Phase 3):
- Full-backbone: 720 sequences / 281.9s = 153.2 seqs/min, mean 9.40s/protein.
- CA-only: 720 sequences / 262.6s = 164.5 seqs/min, mean 8.75s/protein (~7%
  faster, consistent with a smaller per-residue feature tensor).

**Peak memory** (measured with `/usr/bin/time -l`, not captured in the main
batch runs above — added on re-verification; single-protein, 8 seqs/target,
T=0.1, `memory_check/` logs):
| case | length | maximum RSS | peak memory footprint |
|------|--------|-------------|------------------------|
| full-backbone, 3PGK (largest single chain) | 415 | 638.7 MB | 330.8 MB |
| CA-only, 3PGK | 415 | 588.9 MB | 330.1 MB |
| full-backbone, 1TUP (largest total length incl. DNA chains) | 627 | 862.1 MB | 412.3 MB |

Memory scales with sequence length as expected; CA-only is marginally lighter
(fewer backbone atoms per residue) at comparable accuracy cost. All figures
are well within a laptop-class 16GB RAM budget — relevant to the Phase 3
low-resource-inference claim that ProteinMPNN doesn't require a
high-end GPU/server.

**Failure rate**: 0/30 for both models (all structures parsed and completed
inference without error, including the 9 stress-test structures with missing/
gapped residues and the 1TUP DNA complex under explicit chain restriction).

## Files
- `chain_id.jsonl`, `build_chain_id_jsonl.py`
- `fullbackbone/`, `ca_only/` — `seqs/*.fa`, `scores/*.npz`, `terminal_output.log`
- `per_sample_metrics.csv`, `runtime_metrics.csv`, `temperature_sweep_summary.csv`
- `compute_metrics.py`, `make_plots.py`
- `temperature_sweep_comparison.png`, `runtime_comparison.png`
- `completeness_comparison.py` / `completeness_comparison.csv` / `.png` —
  clean vs. stress-test subset comparison (input completeness)
- `memory_check/` — peak memory measurements (`/usr/bin/time -l`)

## Environment
Same as Phase 0/1: conda env `proteinmpnn`, Python 3.10.20, PyTorch 2.12.0,
numpy 2.2.6, Apple M3 / 16GB RAM, CPU only (no CUDA available, and the
official code doesn't use MPS either — see Phase 0 note).

## Deviations from README defaults
- `--sampling_temp "0.1 0.2 0.3"` in a single run (the flag accepts a
  space-separated string and samples `num_seq_per_target` sequences at each
  temperature per the code, `protein_mpnn_run.py:62,324`) rather than three
  separate invocations — functionally identical, fewer model reloads.
- `--num_seq_per_target 8` instead of the example scripts' default of 2, to
  get more stable per-protein recovery/score means.
- `--ca_only` combined with `--jsonl_path` (parsed once in Phase 1) instead of
  `--pdb_path`, which is the only path shown in the README/examples. Confirmed
  via `tied_featurize` (`protein_mpnn_utils.py:197`) that this is a supported,
  equivalent code path — not a hack.

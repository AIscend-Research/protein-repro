# Reproducibility Checklist

Consolidated environment, checkpoint, dataset, and command reference for the
whole reproduction (Phases 0-3). Phase-specific detail lives in each phase's
own README; this is the single-page summary.

## Environment
- OS / hardware: macOS (Darwin), Apple M3, 16 GB RAM
- conda env: `proteinmpnn`, Python 3.10.20
- PyTorch 2.12.0, **CPU only** — the official code's device selection
  (`protein_mpnn_run.py:68`) checks only `torch.cuda.is_available()`, never
  MPS; a naive one-line MPS patch was tested in Phase 3 and crashes on a
  `torch.gather` operation in `gather_nodes`/`gather_edges`
  (`protein_mpnn_utils.py:595-608`) — Apple Silicon GPU acceleration is not
  usable with this codebase as-is (see `phase3_extensions/low_resource/`)
- numpy 2.2.6, scipy 1.15.3, pandas 2.3.3, biopython 1.87, tqdm 4.68.1,
  matplotlib 3.10.9, scikit-learn 1.7.2, torchvision 0.27.0, torchaudio 2.11.0
- Repo git commit: `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57` (unchanged
  throughout the entire reproduction — confirmed via `git rev-parse HEAD`
  before every phase)

## Model checkpoints
| checkpoint | file | params | edges | training noise |
|---|---|---|---|---|
| full-backbone | `vanilla_model_weights/v_48_020.pt` | 1.66M | 48 | 0.20 Å |
| CA-only | `ca_model_weights/v_48_020.pt` | 1.65M | 48 | 0.20 Å |

`v_48_020` used throughout (both models) so the full-backbone vs. CA-only
comparison isn't confounded by a different training noise level. All three
weight sets (vanilla ×4, soluble ×4, ca_only ×3) verified present and
loadable in Phase 0.

## Dataset
30 PDB structures curated from RCSB, `inputs/phase1_dataset/pdbs/` — full
per-structure table (chain IDs, lengths, missing residues, resolution,
clean/stress-test subset) in `reproduction/phase1_data_prep/dataset_summary.csv`.

- 21 clean / 9 stress-test (naturally incomplete or gapped)
- 27 monomers / 3 multimers (1TIM homodimer, 4HHB hetero-tetramer, 1TUP
  hetero-trimer + 2 DNA chains)
- Known caveats (documented in Phase 1 README): 1AKI/1HEL and 1UBI/1UBQ are
  each duplicate protein identities (28 unique sequences, not 30); 1TUP
  requires explicit chain restriction (`assign_fixed_chains.py`) to exclude
  its DNA chains from design.
- Native sequences: `reproduction/phase1_data_prep/native_sequences.fasta`

## Seeds, sampling, commands
- Random seed: **37** throughout (Phase 2 temperature sweep, Phase 3 noise/
  masking sweeps) — chosen arbitrarily once and reused for consistency.
- Sampling temperatures: 0.1 (standard), 0.2, 0.3 — passed as a single
  `--sampling_temp "0.1 0.2 0.3"` argument (samples independently at each).
- Sequences per backbone: 8 per temperature (24 total per protein in the
  3-temperature sweep).
- Exact commands for every experiment are in each phase's README:
  - Phase 0: `reproduction/phase0_smoke_test/README.md`
  - Phase 1: `reproduction/phase1_data_prep/README.md`
  - Phase 2: `reproduction/phase2_reproduction/README.md`
  - Phase 3: `reproduction/phase3_extensions/README.md`

## Headline results
| finding | value |
|---|---|
| Full-backbone recovery (T=0.1, 30 proteins) | 46.4% |
| CA-only recovery (T=0.1, 30 proteins) | 41.3% |
| Recovery vs. temperature | monotonic decrease, both models |
| Recovery vs. Gaussian backbone noise (0 → 0.5 Å) | 46.5% → 23.9% (full-backbone), 40.5% → 25.4% (CA-only); **CA-only more noise-robust past ~0.4 Å** |
| Recovery vs. synthetic residue masking (0 → 30%) | see Phase 3 masking results |
| Buried vs. exposed residue recovery | 42.2% vs. 33.3% (full-backbone) |
| Failure rate (30/30 structures, both models) | 0% |
| Peak memory (largest structure, 1TUP) | 862 MB RSS |
| Peak memory (synthetic 3386-residue stress test) | 2.03 GB footprint, completes in 24.5s, no failure |
| batch_size 1 vs. 8 (num_seq=8, 3PGK) | 641 MB vs. 2258 MB peak RSS (~3.5×) for ~28% speedup |
| CPU-only runtime | ~150-165 sequences/min (Apple M3, 8 cores) |
| ESMFold validation (5 designs, via Kaggle) | 3/5 fold with high confidence & low RMSD (≤1.0 Å) to native despite 22-59% sequence recovery |

## Known limitations
- Small (30-structure, 28-unique-sequence) hand-curated test set, not the
  paper's held-out CATH-clustered test split (not separately downloadable —
  full 16.5 GB training tarball required for `test_clusters.txt`).
- No wet-lab validation performed. Structure-prediction validation (ESMFold)
  was run via the user's Kaggle account
  (`reproduction/phase3_extensions/kaggle_esmfold/`, results in `results/`)
  since this machine has no usable GPU locally — 3 of 5 designed sequences
  fold with high confidence and low CA RMSD to the native structure.
- Buried/exposed analysis uses a CA-CA contact-number proxy, not true DSSP
  solvent accessibility (`mkdssp` not resolvable via conda/Homebrew on this
  osx-arm64 machine).

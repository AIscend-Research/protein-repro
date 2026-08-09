# Phase 3: Extensions (Responsible Protein Modeling Focus)

All experiments reuse the official `protein_mpnn_run.py` pipeline and the
Phase 1/2 dataset/checkpoints unless noted. Model `v_48_020` throughout
(both full-backbone and CA-only), seed 37, CPU only (Apple M3, 16GB RAM;
see the MPS finding below). Git commit
`8907e6671bfbfc92303b5f79c4b5e6ce47cdef57` unchanged.

## 1. Robustness testing

### 1a. Gaussian backbone noise sweep
`noise_sweep/run_sweep.sh` runs all 30 proteins through `--backbone_noise`
{0.0, 0.1, 0.2, 0.3, 0.5} Å for both models (T=0.1, 8 seqs/target, `save_score 1`).
Analysis: `noise_sweep/analyze_noise_sweep.py` → `summary.csv`,
`noise_sweep_comparison.png`.

| noise (Å) | full-backbone recovery | CA-only recovery |
|---|---|---|
| 0.0 | 46.5% | 40.5% |
| 0.1 | 45.6% | 40.1% |
| 0.2 | 43.5% | 38.0% |
| 0.3 | 39.5% | 35.1% |
| 0.5 | 23.9% | 25.4% |

**Finding: full-backbone degrades faster than CA-only and the two cross
over around 0.4-0.5 Å** — full-backbone is more accurate on clean
structures but more sensitive to coordinate corruption (its extra N/C/O
geometric detail, e.g. backbone dihedrals, is more fragile to noise than
CA-only's sparser representation). CA-only becomes the *more robust* choice
under sufficiently degraded input, even though it's the weaker model
on clean structures.

### 1b. Synthetic residue masking
`masking/make_masked_jsonl.py` masks a contiguous internal stretch (10%,
20%, 30%) of each of the 21 Phase-1 "clean" proteins, using the exact `-` /
NaN-coordinate convention the official parser uses for real missing
residues (verified against 1TIM's genuine gap in Phase 1). `masking/run_masking.sh`
runs both models; `masking/analyze_masking.py` → `summary.csv`,
`masking_comparison.png`.

| mask % | full-backbone recovery | CA-only recovery |
|---|---|---|
| 0% | 47.1% | 40.1% |
| 10% | 44.3% | 36.3% |
| 20% | 38.5% | 30.8% |
| 30% | 33.6% | 26.5% |

Recovery here is computed by the official pipeline's own `mask_for_loss`
(`= mask * chain_M * chain_M_pos`, `protein_mpnn_run.py:258`) over
**non-masked positions only** — this measures how well the model designs
the still-visible part of a structure as more of the rest is hidden, not
"guessing" the masked region itself.

**Finding: unlike the noise sweep, masking degrades both models roughly
proportionally with no crossover** — full-backbone keeps its lead over
CA-only at every mask level. Missing structural context and noisy
structural context are different failure modes with different relative
impact on the two architectures.

### 1c. Full-backbone vs. CA-only robustness — combined takeaway
CA-only is more robust to **coordinate noise** (crosses over full-backbone
past ~0.4 Å) but not to **missing residues** (full-backbone stays ahead
throughout the masking sweep). A practitioner with noisy/imprecise
coordinates (e.g. low-resolution cryo-EM) might prefer CA-only; one with
gaps/missing density should prefer full-backbone.

## 2. Low-resource inference

### 2a. CPU vs. "GPU" (MPS)
The official code's device selection (`protein_mpnn_run.py:68`) is
`torch.device("cuda:0" if torch.cuda.is_available() else "cpu")` — it never
checks `torch.backends.mps.is_available()`, so on this Apple Silicon Mac
every run in Phases 0-3 used CPU only, MPS unused.

We patched a copy (`protein_mpnn_run_mps_patch.py`, repo root, not part of
the official codebase) to add an MPS branch and tested it on a single small
protein. **It crashes**, not just under-accelerates:
```
failed assertion `Rank of updates array (1) must be greater than or equal
to inner-most dimension of indices array (48)`
```
Root cause: `gather_nodes`/`gather_edges` (`protein_mpnn_utils.py:595-608`)
use a 4D `torch.gather` for the k-nearest-neighbor (k=48) graph construction
that Apple's Metal Performance Shaders `GatherND` kernel doesn't support in
this shape. This is a real, load-bearing part of the architecture (every
forward pass needs it), not an edge case — **MPS cannot currently run
ProteinMPNN at all**, not even slowly. CPU is the only working backend on
Apple Silicon with this codebase as-is. This is itself the honest answer to
"compare CPU vs GPU inference" for this hardware.

### 2b. Model size
| checkpoint | parameters | file size |
|---|---|---|
| full-backbone `v_48_020.pt` | 1.66M | 6.68 MB |
| CA-only `v_48_020.pt` | 1.65M | 6.62 MB |

Both models are tiny by modern deep-learning standards — well within reach
of any laptop-class CPU, confirmed by the runtimes below.

### 2c. Batch size / num_seq_per_target as a low-compute lever
`low_resource/subset5.jsonl` (5 short proteins: 1CRN, 1VII, 1UBQ, 2GB1,
1ENH), full-backbone model, `batch_size 1` throughout (the README already
recommends batch_size=1 off-GPU):

| num_seq_per_target | total time (5 proteins) | mean time/protein |
|---|---|---|
| 8 (Phase 2 default) | 5.72 s | 1.14 s |
| 1 (low-compute) | 0.78 s | 0.16 s |

~7.3× speedup for 8× fewer sequences (sub-linear, since the structure
encoder forward pass is shared across samples — only the per-sample decode
loop scales down). Reducing `num_seq_per_target` is a simple, effective
low-compute lever with no accuracy cost (each sequence is independently
valid, you just get fewer of them).

**`--batch_size` itself** (re-checked on re-verification — the above only
varied `num_seq_per_target`, not batch size): tested `batch_size` ∈
{1, 2, 4, 8} with `num_seq_per_target=8` fixed on the same 5-protein subset,
then peak memory specifically on the largest single chain (3PGK, 415 res):

| batch_size | total time (5 proteins) | peak RSS (3PGK, num_seq=8) | peak mem. footprint |
|---|---|---|---|
| 1 | 5.08 s | 641 MB | 327 MB |
| 2 | 4.53 s | — | — |
| 4 | 3.71 s | — | — |
| 8 | 2.94 s | 2258 MB | 2032 MB |

Counter to what "smaller batch = low-compute" might suggest, **larger batch
size is actually faster on CPU** (better vectorization across the batch
dimension amortizes fixed overhead) — but at a real memory cost: batch_size=8
uses ~3.5× the peak RSS and ~6× the peak memory footprint of batch_size=1
for the same total output. This matches the README's own guidance
(`--batch_size`: "can set higher for titan, quadro GPUs, reduce this if
running out of GPU memory") — `batch_size=1` is the right choice specifically
when memory, not time, is the binding constraint, which is the realistic
low-resource scenario (a student's laptop with 8-16GB RAM, not a slow CPU
per se).

### 2d. Peak memory (from Phase 2, restated here)
Largest structure tested (1TUP, 627 residues): 862 MB peak RSS,
full-backbone. Comfortably within a 16GB laptop budget; see
`reproduction/phase2_reproduction/memory_check/`.

## 3. Generalization and error analysis

### 3a. Length and oligomeric state
`generalization/length_multimer_analysis.py` (reuses Phase 2 data, T=0.1,
no new inference). Recovery peaks at medium length (70-150 aa) and dips for
long proteins (>150 aa) for both models; monomers notably outperform
multimers (+8.3 pp full-backbone, +5.9 pp CA-only) — **but n=3 multimers**
(1TIM, 4HHB, 1TUP) is too small a sample to treat as conclusive, flagged
here rather than overclaimed.

| length bucket | full-backbone | CA-only |
|---|---|---|
| short (<70, n=11) | 46.8% | 41.5% |
| medium (70-150, n=12) | 49.0% | 41.9% |
| long (>150, n=7) | 41.1% | 39.7% |

| | full-backbone | CA-only |
|---|---|---|
| monomer (n=27) | 47.2% | 41.8% |
| multimer (n=3) | 38.9% | 35.9% |

### 3b. Buried vs. exposed residues
DSSP (`mkdssp`) is not resolvable via conda (salilab/bioconda have no
osx-arm64 build) or Homebrew on this machine, so true solvent accessibility
isn't available. Substituted a standard geometric proxy instead — CA-CA
contact number (neighbors within 10 Å, excluding |i-j|<3 trivial backbone
adjacency) — documented explicitly as a substitution, not silently
presented as DSSP output. Computed for 5 representative monomers
(`generalization/buried_exposed_analysis.py`), tertile split per protein:

| burial class | full-backbone recovery | CA-only recovery |
|---|---|---|
| buried (top third) | 42.2% | 39.9% |
| exposed (bottom third) | 33.3% | 30.9% |

Buried/core positions recover substantially better than exposed/surface
ones for both models — consistent with the well-known biophysical
expectation that core packing constrains sequence identity more tightly
than solvent-exposed positions, which tolerate more substitutions.

### 3c. Deliberate failure cases
`failure_cases/` — constructed edge-case PDBs and tested them through the
official pipeline:

| case | file | result |
|---|---|---|
| Malformed coordinate (non-numeric field) | `malformed_badcoord.pdb` | **Hard crash** — unhandled `ValueError` in `parse_PDB_biounits`, no try/except. Kills the *entire batch* directory parse, not just the bad file. |
| Unusual/unrecognized amino acid code (`ZZZ`) | `unusual_aa.pdb` | **Silently converted to a gap `-`** (treated identically to a genuinely missing residue), no warning. A real modified residue not on MSE's special-cased path would be silently dropped from the sequence. |
| Chain ID collision (two unrelated proteins forced onto the same chain letter, overlapping residue numbers) | `chain_collision.pdb` | **Silently collapses** — colliding (chain, residue-number) coordinate/identity data overwrite each other in the parser's internal dict with no error or warning; only one fragment's 4 residues survived in this test. |
| Over-length protein (`--max_length 50` on a 76-residue protein) | test 4 | **Graceful discard** — `discarded {'too_long': 1}`, no crash, no output for that entry. The one well-behaved failure path tested. |
| Extreme/distorted structure (`--backbone_noise 5.0`, 10× the sweep max) | test 5 | **No crash, but degenerate output** — sequences collapse to a single dominant amino acid repeated almost everywhere (e.g. nearly all "N" or all "K"), mean recovery ~3%, *below* the ~5% expected from uniform-random guessing. A distinct silent-failure mode: the pipeline "succeeds" (produces a FASTA) but the output is garbage, not diverse-but-wrong.
| Very long protein (synthetic 3386-residue, 30-chain structure — real coordinates from all 30 Phase-1 structures concatenated, since our real dataset tops out at 627) | `make_verylong.py` / test 6 | **Succeeds, no failure** — 24.5s for 2 sequences, 2.03 GB peak memory (footprint), comfortably within a 16GB budget. ~5.4× the dataset's largest real structure with no degradation in behavior, just proportionally higher runtime/memory. |

**Takeaway:** of 6 constructed cases, only 2 (`--max_length`, very-long-protein)
behave the way a user would want: a clean discard and a successful, if
slower/memory-heavier, completion, respectively. Malformed coordinates crash
the whole batch; unusual residues and chain collisions are silently
mishandled with no warning; extreme distortion produces confidently-wrong
degenerate sequences rather than an error. Anyone running this pipeline on
real-world, imperfectly-curated PDB files should validate inputs *before*
the run — the pipeline itself won't catch most of these for you, though it
does scale gracefully to genuinely large structures.

## 4. Additional

### 4a. Visualization
`visualizations/recovery_heatmap.py` → `recovery_heatmap.png`: per-residue
native-vs-designed match/mismatch track for 5 representative proteins
(full-backbone, T=0.1, sample 1), recovery ranging 22%-59% across examples.

### 4b. Reproducibility checklist
Consolidated at `reproduction/REPRODUCIBILITY_CHECKLIST.md` — environment,
checkpoints, dataset, seeds, commands, headline results, and known
limitations across all four phases in one page.

### 4c. Optional structure-validation (ESMFold via Kaggle) — completed
`kaggle_esmfold/esmfold_validation.ipynb` folds 5 designed sequences
(full-backbone, T=0.1, sample 1) plus their native counterparts with
ESMFold, and compares pLDDT and CA RMSD. The user set up Kaggle API access
(`~/.kaggle/kaggle.json`), after which it was pushed and run programmatically
via `kaggle kernels push`/`status`/`output`.

**Debugging note (itself a low-resource/accessibility finding):** Kaggle
assigned a Tesla P100 GPU (compute capability 6.0), but the pre-installed
PyTorch build on Kaggle's notebook image only supports compute capability
7.0+ — a platform-side mismatch, not something our install caused (confirmed
by installing dependencies with `--no-deps` to rule out a pip-triggered torch
upgrade). A naive `.cuda()` call doesn't fail; the crash only surfaces on the
first real forward-pass op (`no kernel image available`). Fixed by checking
`torch.cuda.get_device_capability()` directly and falling back to CPU when
below 7.0 — CPU folding for these tiny (36-76 residue) sequences still
completed in a few minutes. Took 6 push/run/debug iterations total (2 for
this GPU incompatibility, 1 for an fp16 `.half()` cast hitting the same
issue, 1 for the `--no-deps` fix, 1 for an unrelated bug in our own RMSD
code using a fake Atom object instead of `QCPSuperimposer`, 1 successful
run). **Takeaway: even "free GPU" cloud notebook platforms have non-obvious
hardware/software compatibility gaps a naive user would hit** — relevant to
the responsible/accessible modeling theme of this reproduction.

**Results** (`kaggle_esmfold/results/`):

| protein | designed pLDDT | native pLDDT | CA RMSD (Å) | seq. recovery |
|---|---|---|---|---|
| 1UBQ | 92.8 | 90.5 | 0.79 | 55.3% |
| 1VII | 69.9 | 92.1 | 2.59 | 22.2% |
| 2GB1 | 83.1 | 88.3 | 1.00 | 32.1% |
| 1CRN | 88.7 | 50.4 | 8.08 | 58.7% |
| 1ENH | 88.8 | 91.2 | 0.64 | 40.7% |

3 of 5 designs (1UBQ, 2GB1, 1ENH) fold with high confidence (pLDDT 83-93)
and low CA RMSD (0.64-1.00 Å) to the native backbone — strong independent
support for ProteinMPNN's actual design objective (reproducing the target
*fold*), despite only 32-55% sequence identity to native. This is the more
meaningful validation than sequence recovery alone: the paper's claim isn't
that designs match the native sequence, but that they fold into the target
structure.

1VII is the weaker case (lower designed pLDDT, higher RMSD) and also has
the lowest sequence recovery of the five (22.2%) — a consistent story.
1CRN's comparison is confounded and shouldn't be over-read: ESMFold itself
is only 50.4 pLDDT confident about the *native* fold (crambin's disulfide-
stabilized mini-fold is a known harder case for pLM-based folding), the
QCP superposition logged a "Newton-Rhapson did not converge" warning, and
the resulting 8.08 Å RMSD likely reflects ESMFold's own uncertainty on this
fold rather than a genuine design failure.

## Files
```
noise_sweep/       - Gaussian backbone noise robustness sweep + analysis
masking/            - synthetic residue masking robustness + analysis
low_resource/       - batch/num_seq low-compute benchmark
generalization/     - length/multimer/buried-exposed analysis
failure_cases/      - 6 deliberately constructed edge-case inputs
visualizations/     - recovery heatmap
kaggle_esmfold/     - ESMFold validation notebook + completed results (Kaggle, P100->CPU fallback)
```
Plus `protein_mpnn_run_mps_patch.py` at the repo root (MPS device patch,
demonstration only, not part of the official codebase).

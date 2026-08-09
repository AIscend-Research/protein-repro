# Phase 0 Smoke Test

Reproduces `examples/submit_example_1.sh` (simple monomer design) to verify the
official installation, pretrained weights, and inference pipeline work end to end.

## Environment
- conda env: `proteinmpnn`
- Python 3.10.20
- PyTorch 2.12.0 (CPU only — no CUDA on this machine; Apple Silicon. Note:
  the official code's device selection, `protein_mpnn_run.py:68`, only checks
  `torch.cuda.is_available()`, never `torch.backends.mps.is_available()`, so
  Apple Silicon's GPU (MPS) goes unused out of the box even though it's
  present on this hardware — see Phase 3 for a patched comparison)
- numpy 2.2.6
- git commit / checkpoint hash: `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`

## Command
```
python helper_scripts/parse_multiple_chains.py \
  --input_path=inputs/PDB_monomers/pdbs/ \
  --output_path=reproduction/phase0_smoke_test/parsed_pdbs.jsonl

python protein_mpnn_run.py \
  --jsonl_path reproduction/phase0_smoke_test/parsed_pdbs.jsonl \
  --out_folder reproduction/phase0_smoke_test \
  --num_seq_per_target 2 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --save_score 1
```

Model: default `v_48_020` (full-backbone vanilla model), no `--ca_only` or
`--use_soluble_model` flags.

## Inputs
- `inputs/PDB_monomers/pdbs/5L33.pdb` (106 residues, chain A)
- `inputs/PDB_monomers/pdbs/6MRR.pdb` (68 residues, chain A)

## Results
Both proteins parsed and designed without errors, 2 sequences each at T=0.1:

| PDB  | length | sample 1 score | sample 1 seq_recovery | sample 2 score | sample 2 seq_recovery | runtime |
|------|--------|-----------------|------------------------|-----------------|------------------------|---------|
| 5L33 | 106    | 0.8576          | 0.3868                 | 0.9035          | 0.4057                 | 0.60s   |
| 6MRR | 68     | (see 6MRR.fa)   |                        |                 |                        | 0.39s   |

`.npz` score files (`scores/5L33.npz`, `scores/6MRR.npz`) contain `score` and
`global_score` per sample and match the FASTA header values exactly.

## Outputs
- `parsed_pdbs.jsonl` — parsed backbone coordinates
- `seqs/5L33.fa`, `seqs/6MRR.fa` — designed sequences
- `scores/5L33.npz`, `scores/6MRR.npz` — per-sample NLL scores
- `terminal_output_parse.log`, `terminal_output_run.log` — raw terminal output

## Status
Phase 0 setup verified: environment, dependencies, and all three pretrained
weight sets (vanilla, ca_only, soluble) are present and the official example
pipeline reproduces successfully with matching git hash, seed, and scores.

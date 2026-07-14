# Phase 1: Data Preparation

## Dataset
30 curated PDB structures downloaded from RCSB (`inputs/phase1_dataset/pdbs/`),
chosen to span monomers, homo-/hetero-oligomers, a range of lengths (36-627
resolved residues), X-ray and NMR structures, and both fully-resolved and
naturally incomplete backbones.

Full stats: `dataset_summary.csv` (built by `analyze_pdbs.py` + `build_summary.py`).
Native sequences for later recovery scoring: `native_sequences.fasta`.

- Total: 30 structures, all parsed and ran through official inference without errors.
- Monomers: 27. Multimers: 3 (1TIM homodimer, 4HHB hetero-tetramer, 1TUP hetero-trimer + DNA).
- **Clean subset (21)**: zero missing/gap residues, protein-only chains.
- **Stress-test subset (9)**: 1PIN, 1R69, 1RIS, 1SHG, 1STN, 1TIM, 1TUP, 2CI2,
  2CRO — each has naturally unresolved/disordered or internally-gapped
  residues. This is a *naturally occurring* stress-test set, distinct from the
  synthetic noise/masking used in Phase 3.
- 3 structures (1BPI, 1FKB, 3CHY) have alternate-location (altloc) atom
  records; both Biopython and the official parser silently resolve to one
  conformer without error (flagged as `has_altloc` in the table, not treated
  as stress-test since the backbone itself is complete).

## Self-audit correction (1TIM reclassified clean -> stress_test)
A later re-check (`verify_formatting.py`) found 1TIM has an internal residue
gap — chain A and B both skip directly from residue 2 to residue 4 — that the
original REMARK-465-based missing-residue count (`analyze_pdbs.py`) reported
as **0**, because 1TIM's PDB file (a 1981 deposition) has no REMARK 465
records at all. REMARK 465 simply isn't populated in enough of the file to
catch it.

The ground-truth signal turned out to be the official parser's own output:
`parse_multiple_chains.py` fills any residue-numbering gap with a `-`
placeholder in `seq_chain_X`, which `tied_featurize` (`protein_mpnn_utils.py:283`)
then converts to `X` before featurization. Scanning `parsed_pdbs.jsonl` for
`-` characters is a more reliable "does ProteinMPNN see a gap here" check than
trusting each file's REMARK 465 completeness, and `build_summary.py` now uses
it (`num_internal_gaps_official_parser` column) as the authoritative signal
alongside REMARK 465. This also explains why 1PIN's REMARK-465 count (10) is
higher than its internal-gap count (5): 5 of its missing residues are
N-terminal truncation (before the resolved numbering starts, so no gap is
created in the parsed sequence) and only 5 are a true internal break.

This does not affect Phase 2's numeric results (recovery/score/perplexity are
computed per protein regardless of subset label) — only this dataset table's
clean/stress-test classification, which is corrected above.

Also verified during this audit: for every one of the 30 structures, the
official parser's total chain-sequence length equals
`biopython_resolved_residues + official_internal_gaps` exactly (checked
programmatically, no exceptions) — including 1PIN, which has a bound
dipeptide ligand (HETATM `ALA`/`PRO` at residues 201-202, far outside the
main chain's 6-163 numbering) that both Biopython and the official parser
correctly exclude from the designed sequence. And 1TUP's 72 REMARK-465
missing residues (23/25/24 across chains A/B/C) are *all* C-terminal
truncation beyond each chain's resolved range (e.g. chain A resolved 94-289,
missing residues start at 290) — consistent with the known disordered
tetramerization/regulatory tail of p53, not internal breaks in the core
domain. This is why 1TUP's `num_internal_gaps_official_parser` (42) comes
entirely from its DNA chains, not its protein chains.

## Dataset composition note: 2 accidental duplicate protein identities
A sequence-identity check across all 30 structures found two pairs with
**identical sequences**: 1AKI and 1HEL are both hen egg-white lysozyme
(129 aa, exact match), and 1UBI and 1UBQ are both human ubiquitin (76 aa,
exact match) — each pair is a separate, legitimate RCSB deposition of the
same well-characterized protein, picked independently when curating a list
of classic small test structures. Both members of each pair parsed and ran
correctly, so this isn't a pipeline error, but it means the 30-structure set
has only **28 unique protein sequences**, slightly under-delivering on the
intended sequence diversity. Left as-is rather than re-curated after the
fact (would invalidate the Phase 2 runs); worth swapping one of each pair for
a different fold if the dataset is extended in Phase 3.

## Tooling
Used the repo's own `helper_scripts/parse_multiple_chains.py` for parsing (no
custom PDB parser written) and `helper_scripts/assign_fixed_chains.py` for
chain selection. Structural stats (missing residues, resolution, chain
composition) were independently computed with Biopython in `analyze_pdbs.py`
since the official parser does not expose that metadata.

## Key failure case found: 1TUP (p53 core domain-DNA complex)
The official parser (`parse_multiple_chains.py`) has no concept of nucleic
acid vs. protein chains — it parsed 1TUP's two DNA chains (E, F, 21 nt each)
as if they were designable protein chains, extending the reported "sequence
length" to 627 instead of the 585 residues across the three protein chains
(A, B, C). Running `protein_mpnn_run.py` on this JSONL with no chain
restriction would silently attempt to redesign DNA as amino acids.

Fix: use `assign_fixed_chains.py --chain_list "A B C"` to explicitly mark
chains A/B/C as designed and E/F as fixed before running inference. Verified
the resulting FASTA header shows `fixed_chains=['E', 'F'], designed_chains=['A', 'B', 'C']`
and a sane `seq_recovery` computed only over the protein chains. See
`spotcheck/1TUP_chain_id.jsonl` and `spotcheck/seqs/1TUP.fa`.

**Takeaway for later phases:** any multi-chain PDB pulled from RCSB must be
checked for non-protein chains (DNA/RNA/ligands parsed as chains) before
batch inference; this is not caught by the official pipeline automatically.

## Chain/position handling spot check (5 proteins)
Verified FASTA headers report the expected `designed_chains`/`fixed_chains`
for: 1UBQ (monomer), 4HHB (hetero-tetramer, all 4 chains designed by default),
1TIM (homodimer, both chains designed), 2CI2 (chain ID is "I", not "A" —
confirms chain-naming isn't assumed; also a stress-test case with 18 missing
residues and still runs cleanly), and 1TUP (DNA chains correctly fixed via
explicit chain selection, see above).

## Files
- `analyze_pdbs.py` / `pdb_stats.json` — independent Biopython structural audit
- `verify_formatting.py` — residue-numbering continuity, altloc, insertion-code check
- `parsed_pdbs.jsonl` — official parser output for all 30 structures
- `build_summary.py` / `dataset_summary.csv` — merged summary table + subset split
- `native_sequences.fasta` — native sequences per structure for recovery scoring
- `terminal_output_parse.log`, `terminal_output_run.log` — raw run logs
- `seqs/*.fa`, `scores/*.npz` — full-batch designed sequences and NLL scores
  (model `v_48_020`, seed 37, T=0.1, 1 seq/target — this is a pipeline
  smoke-run over the whole dataset, not the Phase 2 evaluation itself)
- `spotcheck/` — 1TUP chain-restriction demonstration

## Environment
Same as Phase 0: conda env `proteinmpnn`, Python 3.10.20, PyTorch 2.12.0,
CPU only (no CUDA; the official code never checks MPS — see Phase 0 note),
git commit `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`.

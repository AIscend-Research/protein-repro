"""Build the Phase 1 dataset summary table from pdb_stats.json + parsed_pdbs.jsonl,
and split into clean vs. stress-test subsets."""
import json
import csv

with open("reproduction/phase1_data_prep/pdb_stats.json") as f:
    stats = {r["pdb_id"]: r for r in json.load(f)}

official_chains = {}
official_internal_gaps = {}
with open("reproduction/phase1_data_prep/parsed_pdbs.jsonl") as f:
    for line in f:
        d = json.loads(line)
        official_chains[d["name"]] = d.get("num_of_chains")
        # '-' is the official parser's own internal-gap placeholder (protein_mpnn_utils.py
        # converts it to 'X' at featurization time) -- this is ground truth for what
        # ProteinMPNN sees as a chain break, independent of whether REMARK 465 was
        # populated (older depositions, e.g. 1TIM, often lack REMARK 465 entirely).
        official_internal_gaps[d["name"]] = sum(
            d[k].count("-") for k in d if k.startswith("seq_chain")
        )

# Structures where the official parser picks up non-protein chains (DNA/RNA)
# that are not amino acid chains -- requires explicit chain selection.
NON_PROTEIN_CHAIN_NOTE = {
    "1TUP": "chains E,F are DNA (21 nt each); must restrict design to chains A,B,C",
}

# Has alternate-location (altloc) atom records; Bio.PDB/official parser both
# silently take one conformer. Confirmed via verify_formatting.py.
HAS_ALTLOC = {"1BPI", "1FKB", "3CHY"}

rows = []
for pdb_id, r in sorted(stats.items()):
    is_multimer = r["num_chains"] > 1
    missing_remark465 = r["num_missing_residues"]
    internal_gaps = official_internal_gaps.get(pdb_id, 0)
    note = NON_PROTEIN_CHAIN_NOTE.get(pdb_id, "")
    if pdb_id == "1TIM":
        note = ("2 internal residue gaps (1 per chain, res. 3) not reported in "
                "REMARK 465 -- caught via official parser's '-' gap placeholder, "
                "not REMARK 465 (file has none); legacy 1981 deposition")
    if pdb_id in HAS_ALTLOC:
        note = (note + "; " if note else "") + "has alternate-location (altloc) atoms"
    subset = "stress_test" if (missing_remark465 > 0 or internal_gaps > 0 or pdb_id in NON_PROTEIN_CHAIN_NOTE) else "clean"
    rows.append({
        "pdb_id": pdb_id,
        "chain_ids": ",".join(r["chain_ids"]),
        "num_protein_chains": r["num_chains"],
        "num_chains_official_parser": official_chains.get(pdb_id),
        "total_resolved_residues": r["total_resolved_residues"],
        "num_missing_residues_remark465": missing_remark465,
        "num_internal_gaps_official_parser": internal_gaps,
        "has_altloc": pdb_id in HAS_ALTLOC,
        "resolution_A": r["resolution"],
        "num_nmr_models": r["num_nmr_models"],
        "monomer_or_multimer": "multimer" if is_multimer else "monomer",
        "parses_official_pipeline": True,
        "runs_inference_no_error": True,
        "subset": subset,
        "notes": note,
    })

fieldnames = list(rows[0].keys())
with open("reproduction/phase1_data_prep/dataset_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

n_clean = sum(1 for r in rows if r["subset"] == "clean")
n_stress = sum(1 for r in rows if r["subset"] == "stress_test")
n_multi = sum(1 for r in rows if r["monomer_or_multimer"] == "multimer")
print(f"Total: {len(rows)}  clean: {n_clean}  stress_test: {n_stress}  multimers: {n_multi}")

# Native sequence FASTA (chain-concatenated with '/', matching MPNN's own
# multi-chain FASTA convention) for later sequence-recovery scoring.
with open("reproduction/phase1_data_prep/native_sequences.fasta", "w") as f:
    for pdb_id, r in sorted(stats.items()):
        seq = "/".join(r["chain_seqs"][c] for c in r["chain_ids"])
        f.write(f">{pdb_id} chains={','.join(r['chain_ids'])}\n{seq}\n")

print("Wrote dataset_summary.csv and native_sequences.fasta")

"""Phase 3 generalization: buried vs. exposed residue recovery.

DSSP (mkdssp/dssp binary) is not installed on this machine and isn't
resolvable via conda (salilab/bioconda channels don't have an osx-arm64
build) or Homebrew, so true solvent accessibility isn't available. Substitute
a standard geometric burial proxy instead: CA-CA contact number (count of
other CA atoms within a distance cutoff, excluding close sequence neighbors
to avoid trivial backbone-adjacency inflation) -- higher contact number =
more buried/packed, lower = more surface-exposed. This is a well-established
lightweight substitute when DSSP is unavailable, documented here explicitly
as a substitution rather than silently presented as solvent accessibility.
"""
import json
import csv
import numpy as np

CONTACT_CUTOFF = 10.0  # Angstrom, CA-CA
SEQ_SEP_MIN = 3  # exclude |i-j| < 3 (trivial local backbone contacts)

# Use a handful of representative monomeric clean proteins spanning short/
# medium/long, matching the Phase 1 spot-check style.
EXAMPLE_PROTEINS = ["1UBQ", "2LZM", "1MBN", "3PGK", "2PCY"]

parsed = {}
with open("reproduction/phase1_data_prep/parsed_pdbs.jsonl") as f:
    for line in f:
        d = json.loads(line)
        if d["name"] in EXAMPLE_PROTEINS:
            parsed[d["name"]] = d

native_seqs = {}
with open("reproduction/phase1_data_prep/native_sequences.fasta") as f:
    lines = f.read().strip().split("\n")
for i in range(0, len(lines), 2):
    pdb_id = lines[i][1:].split()[0]
    native_seqs[pdb_id] = lines[i + 1]


def contact_numbers(ca_coords):
    ca = np.array(ca_coords)
    n = len(ca)
    dists = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    contacts = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(n):
            if abs(i - j) >= SEQ_SEP_MIN and dists[i, j] <= CONTACT_CUTOFF:
                contacts[i] += 1
    return contacts


def load_designed_seq(model_folder, pdb_id, sample_idx=0):
    path = f"reproduction/phase2_reproduction/{model_folder}/seqs/{pdb_id}.fa"
    with open(path) as f:
        lines = f.read().strip().split("\n")
    # record 0 = native header+seq, record 1+ = designed samples
    seq = lines[2 * (sample_idx + 1) + 1]
    return seq.replace("/", "")  # single-chain proteins only here


rows = []
for pdb_id in EXAMPLE_PROTEINS:
    d = parsed[pdb_id]
    chain = d["chain_ids"] if "chain_ids" in d else None
    # single-chain proteins in this example set: find the one seq_chain_X key
    chain_letter = [k[-1] for k in d if k.startswith("seq_chain")][0]
    ca_coords = d[f"coords_chain_{chain_letter}"][f"CA_chain_{chain_letter}"]
    contacts = contact_numbers(ca_coords)

    native = native_seqs[pdb_id]
    for model_folder in ["fullbackbone", "ca_only"]:
        designed = load_designed_seq(model_folder, pdb_id, sample_idx=0)
        assert len(designed) == len(native) == len(contacts), \
            f"{pdb_id} {model_folder}: length mismatch {len(designed)} {len(native)} {len(contacts)}"
        for i in range(len(native)):
            rows.append({
                "pdb_id": pdb_id,
                "model": model_folder,
                "position": i,
                "contact_number": int(contacts[i]),
                "match": int(designed[i] == native[i]),
            })

with open("reproduction/phase3_extensions/generalization/buried_exposed_per_residue.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# Tertile split per protein (buried = top third contact number, exposed = bottom third)
from collections import defaultdict
by_protein_contacts = defaultdict(list)
for r in rows:
    if r["model"] == "fullbackbone":  # tertile cutoffs computed once per protein, same for both models
        by_protein_contacts[r["pdb_id"]].append(r["contact_number"])

tertile_cutoffs = {}
for pdb_id, vals in by_protein_contacts.items():
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    tertile_cutoffs[pdb_id] = (vals_sorted[n // 3], vals_sorted[2 * n // 3])

for r in rows:
    lo, hi = tertile_cutoffs[r["pdb_id"]]
    if r["contact_number"] <= lo:
        r["burial_class"] = "exposed"
    elif r["contact_number"] >= hi:
        r["burial_class"] = "buried"
    else:
        r["burial_class"] = "intermediate"

summary = defaultdict(list)
for r in rows:
    if r["burial_class"] in ("buried", "exposed"):
        summary[(r["model"], r["burial_class"])].append(r["match"])

print(f"{'model':<14}{'burial_class':<14}{'n_positions':<14}{'mean_recovery'}")
summary_rows = []
for model in ["fullbackbone", "ca_only"]:
    for cls in ["buried", "exposed"]:
        vals = summary[(model, cls)]
        mean_r = sum(vals) / len(vals)
        summary_rows.append({"model": model, "burial_class": cls,
                              "n_positions": len(vals), "mean_recovery": mean_r})
        print(f"{model:<14}{cls:<14}{len(vals):<14}{mean_r:.4f}")

with open("reproduction/phase3_extensions/generalization/buried_exposed_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
    w.writeheader()
    w.writerows(summary_rows)
print("\nWrote buried_exposed_per_residue.csv and buried_exposed_summary.csv")

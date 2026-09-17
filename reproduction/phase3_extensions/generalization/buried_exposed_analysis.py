"""Phase 3 generalization: buried vs. exposed residue recovery.

Originally used a CA-CA contact-number proxy for burial because the mkdssp
binary has no osx-arm64 conda/Homebrew build. Replaced (on re-verification)
with real per-residue solvent-accessible surface area (SASA) from
`mdtraj.shrake_rupley` (Shrake-Rupley rolling-probe algorithm, the same
algorithm DSSP itself uses internally for accessibility) -- this is a
self-contained Python/C implementation with no external DSSP binary
dependency, so it removes the "proxy, not the real thing" caveat entirely.
Lower SASA = more buried, higher SASA = more solvent-exposed.
"""
import json
import csv
from collections import defaultdict

import mdtraj as md

# Use a handful of representative monomeric clean proteins spanning short/
# medium/long, matching the Phase 1 spot-check style.
EXAMPLE_PROTEINS = ["1UBQ", "2LZM", "1MBN", "3PGK", "2PCY"]
PDB_DIR = "inputs/phase1_dataset/pdbs"

native_seqs = {}
with open("reproduction/phase1_data_prep/native_sequences.fasta") as f:
    lines = f.read().strip().split("\n")
for i in range(0, len(lines), 2):
    pdb_id = lines[i][1:].split()[0]
    native_seqs[pdb_id] = lines[i + 1]


def residue_sasa(pdb_id):
    traj = md.load_pdb(f"{PDB_DIR}/{pdb_id}.pdb")
    prot = traj.atom_slice(traj.topology.select("protein"))
    sasa = md.shrake_rupley(prot, mode="residue")[0]  # nm^2 per residue, single frame
    # mdtraj's "protein" atom selection also picks up non-standard-AA caps
    # (e.g. an ACE acetyl cap on 3PGK) as their own pseudo-residue; drop any
    # residue without a standard one-letter code so indices align 1:1 with
    # the native sequence used elsewhere in the pipeline.
    residues = list(prot.topology.residues)
    keep = [i for i, r in enumerate(residues) if r.code is not None]
    return sasa[keep]


def load_designed_seq(model_folder, pdb_id, sample_idx=0):
    path = f"reproduction/phase2_reproduction/{model_folder}/seqs/{pdb_id}.fa"
    with open(path) as f:
        lines = f.read().strip().split("\n")
    # record 0 = native header+seq, record 1+ = designed samples
    seq = lines[2 * (sample_idx + 1) + 1]
    return seq.replace("/", "")  # single-chain proteins only here


rows = []
for pdb_id in EXAMPLE_PROTEINS:
    sasa = residue_sasa(pdb_id)
    native = native_seqs[pdb_id]
    assert len(sasa) == len(native), \
        f"{pdb_id}: SASA length {len(sasa)} != native length {len(native)}"

    for model_folder in ["fullbackbone", "ca_only"]:
        designed = load_designed_seq(model_folder, pdb_id, sample_idx=0)
        assert len(designed) == len(native), \
            f"{pdb_id} {model_folder}: length mismatch {len(designed)} {len(native)}"
        for i in range(len(native)):
            rows.append({
                "pdb_id": pdb_id,
                "model": model_folder,
                "position": i,
                "sasa_nm2": float(sasa[i]),
                "match": int(designed[i] == native[i]),
            })

with open("reproduction/phase3_extensions/generalization/buried_exposed_per_residue.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# Tertile split per protein (buried = bottom third SASA, exposed = top third SASA)
by_protein_sasa = defaultdict(list)
for r in rows:
    if r["model"] == "fullbackbone":  # tertile cutoffs computed once per protein, same for both models
        by_protein_sasa[r["pdb_id"]].append(r["sasa_nm2"])

tertile_cutoffs = {}
for pdb_id, vals in by_protein_sasa.items():
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    tertile_cutoffs[pdb_id] = (vals_sorted[n // 3], vals_sorted[2 * n // 3])

for r in rows:
    lo, hi = tertile_cutoffs[r["pdb_id"]]
    if r["sasa_nm2"] <= lo:
        r["burial_class"] = "buried"
    elif r["sasa_nm2"] >= hi:
        r["burial_class"] = "exposed"
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

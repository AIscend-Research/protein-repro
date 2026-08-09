"""Phase 1 data prep: scan curated PDB files and record structural stats.

For each PDB: chain IDs, per-chain length (resolved residues), number of
chains, monomer/multimer, number of models (NMR ensembles), resolution,
missing-residue count (from REMARK 465), and native sequence per chain.
Does not use the official ProteinMPNN parser (that is checked separately in
run_mpnn_check.py) -- this is an independent structural audit.
"""
import os
import glob
import json
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

warnings.simplefilter("ignore", PDBConstructionWarning)

PDB_DIR = "inputs/phase1_dataset/pdbs"
OUT_JSON = "reproduction/phase1_data_prep/pdb_stats.json"

STANDARD_AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",  # selenomethionine, common crystallographic substitution
}


def parse_remark465(path):
    """Count missing residues listed in REMARK 465 (disordered/unresolved)."""
    count = 0
    in_table = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("REMARK 465"):
                content = line[10:].strip()
                if content.startswith("M RES C SSSEQI") or content.startswith("RES C SSSEQI"):
                    in_table = True
                    continue
                if in_table and content and not content.startswith("MODELS"):
                    parts = content.split()
                    if len(parts) >= 3:
                        count += 1
    return count


def get_resolution(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("REMARK   2 RESOLUTION."):
                toks = line.split("RESOLUTION.")[1].split()
                for t in toks:
                    try:
                        return float(t)
                    except ValueError:
                        continue
    return None


def get_num_models(path):
    return sum(1 for line in open(path) if line.startswith("MODEL "))


def main():
    parser = PDBParser(QUIET=True)
    records = []
    for path in sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb"))):
        pdb_id = os.path.basename(path).replace(".pdb", "")
        try:
            structure = parser.get_structure(pdb_id, path)
            model = next(iter(structure))
            chains_info = {}
            for chain in model:
                seq = ""
                for res in chain:
                    resname = res.get_resname()
                    if resname in STANDARD_AA_3TO1 and res.id[0] == " ":
                        seq += STANDARD_AA_3TO1[resname]
                if seq:
                    chains_info[chain.id] = {"length": len(seq), "seq": seq}

            num_missing = parse_remark465(path)
            resolution = get_resolution(path)
            num_models = get_num_models(path)

            records.append({
                "pdb_id": pdb_id,
                "num_chains": len(chains_info),
                "chain_ids": list(chains_info.keys()),
                "chain_lengths": {k: v["length"] for k, v in chains_info.items()},
                "chain_seqs": {k: v["seq"] for k, v in chains_info.items()},
                "total_resolved_residues": sum(v["length"] for v in chains_info.values()),
                "num_missing_residues": num_missing,
                "resolution": resolution,
                "num_nmr_models": num_models,
                "is_multimer": len(chains_info) > 1,
                "parse_error": None,
            })
        except Exception as e:
            records.append({
                "pdb_id": pdb_id,
                "parse_error": str(e),
            })

    with open(OUT_JSON, "w") as fh:
        json.dump(records, fh, indent=2)

    print(f"Parsed {len(records)} PDB files -> {OUT_JSON}")
    for r in records:
        if r.get("parse_error"):
            print(f"  {r['pdb_id']}: PARSE ERROR: {r['parse_error']}")
        else:
            print(f"  {r['pdb_id']}: chains={r['num_chains']} lengths={r['chain_lengths']} "
                  f"missing={r['num_missing_residues']} models={r['num_nmr_models']} res={r['resolution']}")


if __name__ == "__main__":
    main()

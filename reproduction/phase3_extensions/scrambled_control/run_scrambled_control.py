"""ESMFold scrambled-sequence control (reviewer-requested addition to Phase 3
kaggle_esmfold validation). Reviewer concern: 3/5 ESMFold "successes" in the
original validation are on textbook demonstration proteins (ubiquitin,
protein G B1, engrailed homeodomain) that are extremely well represented in
ESMFold's own training distribution, so high pLDDT / low RMSD there may
reflect ESMFold's prior over small well-known folds rather than anything
specific to the ProteinMPNN-designed sequence. Missing control: fold a
composition- and length-matched *scrambled* version of each designed
sequence and confirm it does NOT reach comparably high confidence / low
RMSD to the native fold.

Run locally on CPU (no GPU used) -- reuses the native-sequence ESMFold folds
already produced in ../kaggle_esmfold/results/ instead of re-folding them.
"""
import random
import subprocess
import sys

import torch
from transformers import AutoTokenizer, EsmForProteinFolding

KAGGLE_RESULTS = "../kaggle_esmfold/results"

SEQUENCES = {
    "1UBQ": "MTIFVAREDGTTLELEVEPSDTIAELKKKIEEKTGIPPEEQKLIYKGKVLEDEKTLADYNIEEGDTIKLELVPKGG",
    "1VII": "KPTPEEAKKLLGGTVEEFKKLSPEEQEEAYKKNLAK",
    "2GB1": "KKYKVEIEGKNYTGSFTVEAKNKEEAKEKVKEELKKYGVEGEFYFDEENNTFKVKD",
    "1CRN": "TVCCPSKEARDKYLECLKPGTPKEECAKATGCIIIPGTTCPADYPY",
    "1ENH": "APPLTFSAEQRAALDARFARNPELSDEELAALSAELGLPAEQIRAWFAARRAAA",
}

REPORTED = {  # from kaggle_esmfold/README table, designed-sequence results
    "1UBQ": {"plddt": 92.8, "rmsd": 0.79},
    "1VII": {"plddt": 69.9, "rmsd": 2.59},
    "2GB1": {"plddt": 83.1, "rmsd": 1.00},
    "1CRN": {"plddt": 88.7, "rmsd": 8.08},
    "1ENH": {"plddt": 88.8, "rmsd": 0.64},
}

SCRAMBLE_SEED = 37


def scrambled(seq, seed):
    rng = random.Random(seed)
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def main():
    scrambled_seqs = {k: scrambled(v, SCRAMBLE_SEED) for k, v in SEQUENCES.items()}
    for k in SEQUENCES:
        assert sorted(scrambled_seqs[k]) == sorted(SEQUENCES[k])  # same composition
        assert scrambled_seqs[k] != SEQUENCES[k]
    print("Scrambled sequences (composition- and length-matched):")
    for k, v in scrambled_seqs.items():
        print(f"  {k}: {v}")

    print("Loading ESMFold (facebook/esmfold_v1) on CPU -- this downloads "
          "~2.7GB on first run and can take a few minutes to load...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model.trunk.set_chunk_size(64)
    model.eval()
    print("Model loaded.")

    def fold(seq):
        tok = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        with torch.no_grad():
            out = model(tok["input_ids"])
        mean_plddt = out["plddt"][0, :, 1].mean().item() * 100
        pdb_str = model.output_to_pdb(out)[0]
        return mean_plddt, pdb_str

    results = {}
    for pdb_id, seq in scrambled_seqs.items():
        print(f"Folding {pdb_id} scrambled ({len(seq)} residues)...")
        plddt, pdb_str = fold(seq)
        out_path = f"{pdb_id}_scrambled_esmfold.pdb"
        with open(out_path, "w") as f:
            f.write(pdb_str)
        results[pdb_id] = {"plddt_scrambled": plddt}
        print(f"  {pdb_id} scrambled pLDDT = {plddt:.1f}")

    # CA RMSD of scrambled-fold vs. the already-computed native-sequence fold
    from Bio.PDB import PDBParser
    from Bio.PDB.qcprot import QCPSuperimposer
    import numpy as np

    parser = PDBParser(QUIET=True)

    def get_ca_coords(path):
        s = parser.get_structure("x", path)
        return np.array([res["CA"].coord for res in s[0]["A"] if "CA" in res])

    print(f"\n{'pdb_id':<8}{'plddt_scrambled':<17}{'plddt_designed':<16}{'rmsd_scrambled_A':<18}{'rmsd_designed_A'}")
    rows = []
    for pdb_id in SEQUENCES:
        ca_scrambled = get_ca_coords(f"{pdb_id}_scrambled_esmfold.pdb")
        ca_native = get_ca_coords(f"{KAGGLE_RESULTS}/{pdb_id}_native_esmfold.pdb")
        n = min(len(ca_native), len(ca_scrambled))
        sup = QCPSuperimposer()
        sup.set(ca_native[:n], ca_scrambled[:n])
        sup.run()
        rmsd = sup.get_rms()
        rep = REPORTED[pdb_id]
        print(f"{pdb_id:<8}{results[pdb_id]['plddt_scrambled']:<17.1f}{rep['plddt']:<16.1f}"
              f"{rmsd:<18.2f}{rep['rmsd']:.2f}")
        rows.append({
            "pdb_id": pdb_id, "plddt_scrambled": results[pdb_id]["plddt_scrambled"],
            "plddt_designed_reported": rep["plddt"], "rmsd_scrambled_to_native_fold_A": rmsd,
            "rmsd_designed_to_native_fold_A_reported": rep["rmsd"],
        })

    import csv
    with open("scrambled_control_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nWrote scrambled_control_results.csv")


if __name__ == "__main__":
    main()

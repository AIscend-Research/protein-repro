"""Phase 3 robustness: synthetically mask a contiguous internal stretch of
each clean protein's chains to simulate incomplete experimental structures,
using the exact same '-' / NaN-coordinate convention the official parser
uses for naturally missing residues (verified against 1TIM's real gap in
Phase 1). Only the 21 Phase-1 'clean' proteins are used as the base, so the
masked fraction is precisely controlled (no pre-existing natural gaps).
"""
import json
import random

MASK_FRACTIONS = [0.10, 0.20, 0.30]
MARGIN = 5  # keep the masked stretch away from chain termini by this many residues
SEED = 42

clean_ids = set()
with open("reproduction/phase1_data_prep/dataset_summary.csv") as f:
    import csv
    for r in csv.DictReader(f):
        if r["subset"] == "clean":
            clean_ids.add(r["pdb_id"])

entries = {}
with open("reproduction/phase1_data_prep/parsed_pdbs.jsonl") as f:
    for line in f:
        d = json.loads(line)
        if d["name"] in clean_ids:
            entries[d["name"]] = d

print(f"Base: {len(entries)} clean proteins")

for frac in MASK_FRACTIONS:
    rng = random.Random(SEED)
    out_path = f"reproduction/phase3_extensions/masking/masked_pdbs_frac{int(frac*100)}.jsonl"
    n_masked_total = 0
    n_total = 0
    with open(out_path, "w") as out:
        for name in sorted(entries):
            d = json.loads(json.dumps(entries[name]))  # deep copy
            chain_keys = [k[-1] for k in d if k.startswith("seq_chain")]
            for c in chain_keys:
                seq = d[f"seq_chain_{c}"]
                L = len(seq)
                mask_len = max(1, round(frac * L))
                lo = MARGIN
                hi = L - mask_len - MARGIN
                if hi <= lo:
                    start = max(0, (L - mask_len) // 2)
                else:
                    start = rng.randint(lo, hi)
                end = start + mask_len

                new_seq = seq[:start] + "-" * mask_len + seq[end:]
                d[f"seq_chain_{c}"] = new_seq
                n_masked_total += mask_len
                n_total += L

                for atom_key, coord_list in d[f"coords_chain_{c}"].items():
                    for i in range(start, end):
                        coord_list[i] = [float("nan"), float("nan"), float("nan")]

            # rebuild the top-level concatenated 'seq' field consistently
            d["seq"] = "".join(d[f"seq_chain_{c}"] for c in sorted(chain_keys))
            out.write(json.dumps(d) + "\n")
    print(f"frac={frac}: wrote {out_path}, "
          f"{n_masked_total}/{n_total} residues masked ({n_masked_total/n_total:.1%})")

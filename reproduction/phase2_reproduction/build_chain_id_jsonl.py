"""Build a chain_id_jsonl designating all protein chains as designed for every
structure in the Phase 1 dataset, except 1TUP where DNA chains E/F are fixed
(see reproduction/phase1_data_prep/README.md)."""
import json

SPECIAL_FIXED = {"1TUP": ["E", "F"]}

out = {}
with open("reproduction/phase1_data_prep/parsed_pdbs.jsonl") as f:
    for line in f:
        d = json.loads(line)
        name = d["name"]
        all_chains = [k[-1] for k in d if k.startswith("seq_chain")]
        fixed = SPECIAL_FIXED.get(name, [])
        designed = [c for c in all_chains if c not in fixed]
        out[name] = [designed, fixed]

with open("reproduction/phase2_reproduction/chain_id.jsonl", "w") as f:
    f.write(json.dumps(out) + "\n")

print(f"Wrote chain assignments for {len(out)} structures")
print("1TUP:", out["1TUP"])

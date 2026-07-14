"""Phase 3 failure-case: construct a genuinely very-long synthetic 'protein'
by concatenating the first chain of all 30 real Phase-1 structures into one
multi-chain design entry (real coordinate data, no fabrication -- just many
unrelated real chains combined into one JSONL record to stress-test scaling
behavior). Distinct from the --max_length test, which only exercises the
discard path on an ordinary-length protein."""
import json
import string

ALPHABET = list(string.ascii_uppercase + string.ascii_lowercase)  # 52 available

entries = []
with open("reproduction/phase1_data_prep/parsed_pdbs.jsonl") as f:
    for line in f:
        entries.append(json.loads(line))

merged = {"name": "SYNTHETIC_VERYLONG", "num_of_chains": 0}
letter_idx = 0
total_len = 0
for d in entries:
    first_chain_key = sorted(k for k in d if k.startswith("seq_chain"))[0]
    orig_letter = first_chain_key[-1]
    new_letter = ALPHABET[letter_idx]
    letter_idx += 1
    merged[f"seq_chain_{new_letter}"] = d[f"seq_chain_{orig_letter}"]
    merged[f"coords_chain_{new_letter}"] = {
        f"{atom}_chain_{new_letter}": coords
        for key, coords in d[f"coords_chain_{orig_letter}"].items()
        for atom in [key.split("_chain_")[0]]
    }
    total_len += len(d[f"seq_chain_{orig_letter}"])
    merged["num_of_chains"] += 1

merged["seq"] = "".join(merged[f"seq_chain_{ALPHABET[i]}"] for i in range(letter_idx))

with open("reproduction/phase3_extensions/failure_cases/verylong_synthetic.jsonl", "w") as f:
    f.write(json.dumps(merged) + "\n")

print(f"Built synthetic mega-structure: {letter_idx} chains, {total_len} total residues")
print(f"(for comparison: largest real structure in the dataset is 1TUP at 627 residues total)")

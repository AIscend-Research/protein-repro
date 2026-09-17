#!/bin/bash
# Phase 3, deliberate failure cases. Six constructed edge-case inputs pushed
# through the official pipeline to characterise how it fails (crash / silent
# corruption / graceful discard / degenerate-but-successful output).
#
#   1  malformed coordinate field   -> expected: hard crash in the parser
#   2  unusual residue code (ZZZ)   -> expected: silently becomes a gap "-"
#   3  chain-ID collision           -> expected: silent overwrite, residues lost
#   4  --max_length 50 on 1UBQ (76) -> expected: graceful discard
#   5  --backbone_noise 5.0         -> expected: runs, degenerate output
#   6  synthetic 3386-residue input  -> expected: succeeds, slower/heavier
#
# NOTE: tests 1-3 are *expected* to fail or misbehave. `set -e` is deliberately
# NOT used, and each step is guarded so the script runs to completion; the exit
# codes are the result being measured. Everything is logged for the analyzer.
#
# Usage: bash reproduction/phase3_extensions/failure_cases/run_failure_cases.sh
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
# Prefer the pinned conda env; fall back to whatever python is active so the
# script still runs in a venv/CI shell without conda on PATH.
PYBIN=$(conda run -n proteinmpnn which python 2>/dev/null || command -v python3)
FC=reproduction/phase3_extensions/failure_cases
mkdir -p "$FC/_tmp"

parse_one () {  # pdb_basename, output_jsonl, log
  local pdb=$1 out=$2 log=$3
  local dir="$FC/_tmp/$(basename "$pdb" .pdb)"
  rm -rf "$dir" && mkdir -p "$dir"
  cp "$FC/pdbs/$pdb" "$dir/"
  $PYBIN helper_scripts/parse_multiple_chains.py \
    --input_path "$dir" --output_path "$out" > "$log" 2>&1
  echo "exit_code=$?" >> "$log"
}

echo "=== test 1: malformed coordinate field (expect parser crash) ==="
parse_one malformed_badcoord.pdb "$FC/test1_parsed.jsonl" "$FC/test1_malformed.log"
tail -3 "$FC/test1_malformed.log"

echo "=== test 2: unusual amino acid code ZZZ (expect silent gap) ==="
parse_one unusual_aa.pdb "$FC/test2_parsed.jsonl" "$FC/test2_unusual_aa.log"
tail -3 "$FC/test2_unusual_aa.log"

echo "=== test 3: chain ID collision (expect silent overwrite) ==="
parse_one chain_collision.pdb "$FC/test3_parsed.jsonl" "$FC/test3_chain_collision.log"
tail -3 "$FC/test3_chain_collision.log"

# --- single-protein JSONL for tests 4 and 5 ---------------------------------
$PYBIN - <<'PY'
import json
src = "reproduction/phase1_data_prep/parsed_pdbs.jsonl"
dst = "reproduction/phase3_extensions/failure_cases/1UBQ_parsed.jsonl"
with open(src) as f, open(dst, "w") as out:
    n = sum(out.write(l) is not None for l in f if json.loads(l)["name"] == "1UBQ")
assert n == 1, f"expected 1 entry for 1UBQ, got {n}"
print(f"wrote {dst}")
PY

echo "=== test 4: --max_length 50 on a 76-residue protein (expect discard) ==="
$PYBIN protein_mpnn_run.py \
  --jsonl_path "$FC/1UBQ_parsed.jsonl" \
  --out_folder "$FC/test4_maxlength" \
  --model_name v_48_020 \
  --num_seq_per_target 2 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --max_length 50 > "$FC/test4_maxlength.log" 2>&1
echo "exit_code=$?" >> "$FC/test4_maxlength.log"
grep -i "discard" "$FC/test4_maxlength.log" || tail -3 "$FC/test4_maxlength.log"

echo "=== test 5: --backbone_noise 5.0, 10x the sweep max (expect degenerate) ==="
$PYBIN protein_mpnn_run.py \
  --jsonl_path "$FC/1UBQ_parsed.jsonl" \
  --out_folder "$FC/test5_extremenoise" \
  --model_name v_48_020 \
  --num_seq_per_target 8 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --backbone_noise 5.0 > "$FC/test5_extremenoise.log" 2>&1
echo "exit_code=$?" >> "$FC/test5_extremenoise.log"
tail -3 "$FC/test5_extremenoise.log"

echo "=== test 6: synthetic 3386-residue / 30-chain structure (expect success) ==="
$PYBIN "$FC/make_verylong.py"
# All 30 synthetic chains are designed; the chain_id JSONL lists them explicitly
# rather than relying on the default, so the record is unambiguous.
$PYBIN - <<'PY'
import json, string
FC = "reproduction/phase3_extensions/failure_cases"
with open(f"{FC}/verylong_synthetic.jsonl") as f:
    entry = json.loads(f.readline())
alphabet = list(string.ascii_uppercase + string.ascii_lowercase)
chains = [c for c in alphabet if f"seq_chain_{c}" in entry]
with open(f"{FC}/verylong_chain_id.jsonl", "w") as f:
    f.write(json.dumps({entry["name"]: [chains, []]}) + "\n")
print(f"wrote verylong_chain_id.jsonl: {len(chains)} designed chains")
PY
/usr/bin/time -l $PYBIN protein_mpnn_run.py \
  --jsonl_path "$FC/verylong_synthetic.jsonl" \
  --chain_id_jsonl "$FC/verylong_chain_id.jsonl" \
  --out_folder "$FC/test6_verylong" \
  --model_name v_48_020 \
  --num_seq_per_target 2 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --max_length 40000 > "$FC/test6_verylong.log" 2>&1
echo "exit_code=$?" >> "$FC/test6_verylong.log"
grep -E "generated in|maximum resident set size|peak memory footprint" "$FC/test6_verylong.log" || true

rm -rf "$FC/_tmp"
echo "ALL_DONE — now run: python reproduction/phase3_extensions/failure_cases/analyze_failure_cases.py"

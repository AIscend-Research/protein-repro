#!/bin/bash
# Phase 2 peak-memory measurement. Runs single-protein inference under
# `/usr/bin/time -l` (macOS) to capture maximum RSS and peak memory footprint,
# which the batch runs in the Phase 2 README do not report.
#
# Cases: 3PGK (415 res, largest single chain) full-backbone + CA-only, and
# 1TUP (627 res incl. DNA chains) full-backbone.
#
# Usage: bash reproduction/phase2_reproduction/memory_check/run_memory_check.sh
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
# Prefer the pinned conda env; fall back to whatever python is active so the
# script still runs in a venv/CI shell without conda on PATH.
PYBIN=$(conda run -n proteinmpnn which python 2>/dev/null || command -v python3)
OUT=reproduction/phase2_reproduction/memory_check

# Single-protein JSONLs, filtered out of the Phase 1 parse so no re-parsing
# (and no chance of a different parser invocation) is involved.
for pdb in 3PGK 1TUP; do
  $PYBIN - "$pdb" <<'PY'
import json, sys
pdb = sys.argv[1]
src = "reproduction/phase1_data_prep/parsed_pdbs.jsonl"
dst = f"reproduction/phase2_reproduction/memory_check/{pdb}_parsed.jsonl"
with open(src) as f, open(dst, "w") as out:
    n = 0
    for line in f:
        if json.loads(line)["name"] == pdb:
            out.write(line)
            n += 1
assert n == 1, f"expected exactly 1 entry for {pdb}, got {n}"
print(f"wrote {dst}")
PY
done

run_case () {  # name, jsonl, out_folder, extra args...
  local name=$1 jsonl=$2 out_folder=$3
  shift 3
  echo "=== memory check: $name ==="
  /usr/bin/time -l $PYBIN protein_mpnn_run.py \
    --jsonl_path "$jsonl" \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder "$out_folder" \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    "$@" > "$OUT/${name}.log" 2>&1
  grep -E "maximum resident set size|peak memory footprint" "$OUT/${name}.log" || true
}

# Output folder names match the directories already committed from the original
# run, so a re-run lands on top of the same paths.
run_case fullbackbone_3PGK "$OUT/3PGK_parsed.jsonl" "$OUT/fullbackbone_3PGK_seqs"
run_case ca_only_3PGK      "$OUT/3PGK_parsed.jsonl" "$OUT/ca_only" --ca_only
run_case fullbackbone_1TUP "$OUT/1TUP_parsed.jsonl" "$OUT/fullbackbone_1TUP_seqs"

echo "ALL_DONE — parse with reproduction/make_paper_tables.py"

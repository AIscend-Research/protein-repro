#!/bin/bash
# Phase 3, low-resource inference benchmark.
#
#   1. num_seq_per_target 8 (Phase 2 default) vs. 1 (low-compute)  -> wall time
#   2. batch_size {1,2,4,8} at num_seq_per_target 8                -> wall time
#   3. peak memory at batch_size 1 vs. 8 on the largest single chain (3PGK)
#
# All on the 5-short-protein subset (1CRN, 1VII, 1UBQ, 2GB1, 1ENH), full-backbone
# v_48_020, T=0.1, seed 37. Logs are what analyze_low_resource.py reads.
#
# Usage: bash reproduction/phase3_extensions/low_resource/run_low_resource.sh
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"
# Prefer the pinned conda env; fall back to whatever python is active so the
# script still runs in a venv/CI shell without conda on PATH.
PYBIN=$(conda run -n proteinmpnn which python 2>/dev/null || command -v python3)
OUT=reproduction/phase3_extensions/low_resource
SUBSET="1CRN 1VII 1UBQ 2GB1 1ENH"

# --- subset JSONL (filtered out of the Phase 1 parse, no re-parsing) ---------
$PYBIN - "$SUBSET" <<'PY'
import json, sys
keep = set(sys.argv[1].split())
src = "reproduction/phase1_data_prep/parsed_pdbs.jsonl"
dst = "reproduction/phase3_extensions/low_resource/subset5.jsonl"
with open(src) as f, open(dst, "w") as out:
    names = []
    for line in f:
        name = json.loads(line)["name"]
        if name in keep:
            out.write(line)
            names.append(name)
missing = keep - set(names)
assert not missing, f"missing from parsed_pdbs.jsonl: {sorted(missing)}"
print(f"wrote {dst}: {names}")
PY

run_timed () {  # label, out_folder, num_seq, batch_size
  local label=$1 out_folder=$2 num_seq=$3 batch=$4
  echo "=== $label (num_seq=$num_seq, batch_size=$batch) ==="
  # `time` on the whole invocation: the per-protein "generated in Ns" lines in
  # the log exclude model load, which is exactly the fixed overhead a
  # low-resource user pays, so both are recorded.
  { time $PYBIN protein_mpnn_run.py \
      --jsonl_path "$OUT/subset5.jsonl" \
      --out_folder "$out_folder" \
      --model_name v_48_020 \
      --num_seq_per_target "$num_seq" \
      --sampling_temp "0.1" \
      --seed 37 \
      --batch_size "$batch" ; } > "$OUT/${label}.log" 2>&1
  grep -E "real|generated in" "$OUT/${label}.log" | tail -3 || true
}

# --- 1. num_seq_per_target lever --------------------------------------------
run_timed highcompute_nseq8 "$OUT/highcompute_nseq8" 8 1
run_timed lowcompute_nseq1  "$OUT/lowcompute_nseq1"  1 1

# --- 2. batch_size lever -----------------------------------------------------
for bs in 1 2 4 8; do
  run_timed "batchsize${bs}" "$OUT/batchsize${bs}" 8 "$bs"
done

# --- 3. peak memory, batch_size 1 vs. 8, on the largest single chain ---------
# 3PGK (415 res) rather than the short subset: memory is length-driven, and this
# is the same structure the Phase 2 memory check uses, so the numbers compare.
MEMJSONL=reproduction/phase2_reproduction/memory_check/3PGK_parsed.jsonl
if [ ! -f "$MEMJSONL" ]; then
  echo "!! $MEMJSONL missing — run reproduction/phase2_reproduction/memory_check/run_memory_check.sh first"
  exit 1
fi
for bs in 1 8; do
  echo "=== mem_bs${bs} (3PGK, num_seq=8, batch_size=$bs) ==="
  /usr/bin/time -l $PYBIN protein_mpnn_run.py \
    --jsonl_path "$MEMJSONL" \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder "$OUT/mem_bs${bs}" \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size "$bs" > "$OUT/mem_bs${bs}.log" 2>&1
  grep -E "maximum resident set size|peak memory footprint" "$OUT/mem_bs${bs}.log" || true
done

echo "ALL_DONE — now run: python reproduction/phase3_extensions/low_resource/analyze_low_resource.py"

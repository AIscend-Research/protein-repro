#!/bin/bash
# End-to-end reproduction driver: Phase 0 (smoke test) through Phase 3
# (extensions) plus the paper table dump, in dependency order, from a clean
# checkout. Every step shells out to the same scripts documented in the
# per-phase READMEs — this file adds ordering, not new behaviour.
#
#   bash reproduction/run_all.sh            # everything (~30-45 min, CPU)
#   bash reproduction/run_all.sh 2 3        # only the named phases
#
# Prerequisites (see REPRODUCIBILITY_CHECKLIST.md):
#   - conda env `proteinmpnn` created and populated
#   - inputs/phase1_dataset/pdbs/ populated with the 30 curated PDBs
#   - model weights present (vanilla_model_weights/, ca_model_weights/)
#
# Phase 3's ESMFold structure validation is NOT run here: it needs a GPU and a
# Kaggle account. See phase3_extensions/kaggle_esmfold/README notes.
set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
# Prefer the pinned conda env; fall back to whatever python is active so the
# script still runs in a venv/CI shell without conda on PATH.
PYBIN=$(conda run -n proteinmpnn which python 2>/dev/null || command -v python3)
PHASES="${*:-0 1 2 3 tables}"

echo "repo:   $REPO_ROOT"
echo "commit: $(git rev-parse HEAD)"
echo "python: $PYBIN"
echo "phases: $PHASES"

has_phase () { [[ " $PHASES " == *" $1 "* ]]; }
banner () { echo; echo "############ $* ############"; }

# --- Phase 0: smoke test -----------------------------------------------------
if has_phase 0; then
  banner "PHASE 0 — smoke test on the repo's own example monomers"
  $PYBIN helper_scripts/parse_multiple_chains.py \
    --input_path=inputs/PDB_monomers/pdbs/ \
    --output_path=reproduction/phase0_smoke_test/parsed_pdbs.jsonl
  $PYBIN protein_mpnn_run.py \
    --jsonl_path reproduction/phase0_smoke_test/parsed_pdbs.jsonl \
    --out_folder reproduction/phase0_smoke_test \
    --num_seq_per_target 2 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    --save_score 1 2>&1 | tee reproduction/phase0_smoke_test/terminal_output.log
fi

# --- Phase 1: data preparation ----------------------------------------------
if has_phase 1; then
  banner "PHASE 1 — parse and characterise the 30-structure dataset"
  $PYBIN helper_scripts/parse_multiple_chains.py \
    --input_path=inputs/phase1_dataset/pdbs/ \
    --output_path=reproduction/phase1_data_prep/parsed_pdbs.jsonl
  $PYBIN reproduction/phase1_data_prep/analyze_pdbs.py
  $PYBIN reproduction/phase1_data_prep/build_summary.py
  $PYBIN reproduction/phase1_data_prep/verify_formatting.py
fi

# --- Phase 2: main reproduction ---------------------------------------------
if has_phase 2; then
  banner "PHASE 2 — full-backbone vs. CA-only, temperature sweep"
  $PYBIN reproduction/phase2_reproduction/build_chain_id_jsonl.py

  for model in fullbackbone ca_only; do
    CA_FLAG=""
    [ "$model" = "ca_only" ] && CA_FLAG="--ca_only"
    echo "=== $model ==="
    $PYBIN protein_mpnn_run.py $CA_FLAG \
      --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
      --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
      --out_folder reproduction/phase2_reproduction/$model \
      --model_name v_48_020 \
      --num_seq_per_target 8 \
      --sampling_temp "0.1 0.2 0.3" \
      --seed 37 \
      --batch_size 1 \
      --save_score 1 \
      > reproduction/phase2_reproduction/$model/terminal_output.log 2>&1
  done

  $PYBIN reproduction/phase2_reproduction/compute_metrics.py
  $PYBIN reproduction/phase2_reproduction/completeness_comparison.py
  $PYBIN reproduction/phase2_reproduction/make_plots.py
  bash reproduction/phase2_reproduction/memory_check/run_memory_check.sh
fi

# --- Phase 3: extensions -----------------------------------------------------
if has_phase 3; then
  banner "PHASE 3a — Gaussian backbone noise sweep"
  bash reproduction/phase3_extensions/noise_sweep/run_sweep.sh
  $PYBIN reproduction/phase3_extensions/noise_sweep/analyze_noise_sweep.py

  banner "PHASE 3b — synthetic residue masking"
  $PYBIN reproduction/phase3_extensions/masking/make_masked_jsonl.py
  bash reproduction/phase3_extensions/masking/run_masking.sh
  $PYBIN reproduction/phase3_extensions/masking/analyze_masking.py

  banner "PHASE 3c — low-resource inference benchmark"
  bash reproduction/phase3_extensions/low_resource/run_low_resource.sh
  $PYBIN reproduction/phase3_extensions/low_resource/analyze_low_resource.py

  banner "PHASE 3d — generalization / error analysis"
  $PYBIN reproduction/phase3_extensions/generalization/length_multimer_analysis.py
  $PYBIN reproduction/phase3_extensions/generalization/buried_exposed_analysis.py

  banner "PHASE 3e — deliberate failure cases (some steps are expected to fail)"
  bash reproduction/phase3_extensions/failure_cases/run_failure_cases.sh
  $PYBIN reproduction/phase3_extensions/failure_cases/analyze_failure_cases.py

  banner "PHASE 3f — visualizations"
  $PYBIN reproduction/phase3_extensions/visualizations/recovery_heatmap.py
fi

# --- Paper tables ------------------------------------------------------------
if has_phase tables; then
  banner "TABLES — regenerate every paper table from the result CSVs"
  $PYBIN reproduction/make_paper_tables.py
fi

banner "ALL_DONE"

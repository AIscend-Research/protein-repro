#!/bin/bash
# Multi-seed re-run of the noise levels around the reported full-backbone /
# CA-only crossover (~0.4-0.5 A), to check whether the crossover survives
# noise-realization variance. Seed 37 at noise 0.5 already exists in
# ../noise_sweep/{fullbackbone,ca_only}_noise0.5 and is reused, not rerun.
set -e
cd /Users/christopherhuang/Documents/GitHub/mlrc-ProteinMPNN
PYBIN=$(conda run -n proteinmpnn which python)
OUTBASE=reproduction/phase3_extensions/robustness_addenda/seed_sweep
mkdir -p $OUTBASE
NOISE_LEVELS="0.4 0.45 0.5 0.6"
SEEDS="37 38 39"

for noise in $NOISE_LEVELS; do
  for seed in $SEEDS; do
    if [ "$noise" = "0.5" ] && [ "$seed" = "37" ]; then
      echo "=== skip noise=0.5 seed=37 (reuse existing noise_sweep run) ==="
      continue
    fi
    for model_flag in "fullbackbone:" "ca_only:--ca_only"; do
      name="${model_flag%%:*}"
      flag="${model_flag##*:}"
      outdir="${OUTBASE}/${name}_noise${noise}_seed${seed}"
      if [ -f "${outdir}/seqs/1AKI.fa" ]; then
        echo "=== skip ${name} noise=${noise} seed=${seed} (already done) ==="
        continue
      fi
      echo "=== ${name}, noise=${noise}, seed=${seed} ==="
      $PYBIN protein_mpnn_run.py $flag \
        --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
        --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
        --out_folder "$outdir" \
        --model_name v_48_020 \
        --num_seq_per_target 8 \
        --sampling_temp "0.1" \
        --seed $seed \
        --batch_size 1 \
        --backbone_noise $noise \
        --save_score 1 > "${outdir}.log" 2>&1
    done
  done
done
echo "SEED_SWEEP_ALL_DONE"

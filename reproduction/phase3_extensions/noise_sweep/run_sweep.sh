#!/bin/bash
set -e
cd /Users/christopherhuang/Documents/GitHub/ProteinMPNN
PYBIN=$(conda run -n proteinmpnn which python)
NOISE_LEVELS="0.0 0.1 0.2 0.3 0.5"

for noise in $NOISE_LEVELS; do
  echo "=== full-backbone, noise=$noise ==="
  $PYBIN protein_mpnn_run.py \
    --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder reproduction/phase3_extensions/noise_sweep/fullbackbone_noise${noise} \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    --backbone_noise $noise \
    --save_score 1 > reproduction/phase3_extensions/noise_sweep/fullbackbone_noise${noise}.log 2>&1

  echo "=== ca_only, noise=$noise ==="
  $PYBIN protein_mpnn_run.py --ca_only \
    --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder reproduction/phase3_extensions/noise_sweep/ca_only_noise${noise} \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    --backbone_noise $noise \
    --save_score 1 > reproduction/phase3_extensions/noise_sweep/ca_only_noise${noise}.log 2>&1
done
echo "ALL_DONE"

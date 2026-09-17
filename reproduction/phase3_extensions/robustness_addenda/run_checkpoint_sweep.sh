#!/bin/bash
# Tests whether the noise-crossover point is a property of the CA-only
# architecture, or of the gap between the checkpoint's training noise and
# the sweep's test noise. CA-only ships v_48_002/010/020 (trained with
# 0.02/0.10/0.20 A Gaussian noise); the main sweep only used v_48_020.
# Re-runs the coordinate-noise sweep with the other two CA-only checkpoints,
# holding the full-backbone v_48_020 curve fixed (already in
# ../noise_sweep/fullbackbone_noise*) as the comparison baseline.
set -e
cd /Users/christopherhuang/Documents/GitHub/mlrc-ProteinMPNN
PYBIN=$(conda run -n proteinmpnn which python)
OUTBASE=reproduction/phase3_extensions/robustness_addenda/checkpoint_sweep
mkdir -p $OUTBASE
NOISE_LEVELS="0.0 0.2 0.3 0.4 0.5"
CKPTS="v_48_002 v_48_010"

for ckpt in $CKPTS; do
  for noise in $NOISE_LEVELS; do
    outdir="${OUTBASE}/ca_only_${ckpt}_noise${noise}"
    if [ -f "${outdir}/seqs/1AKI.fa" ]; then
      echo "=== skip ${ckpt} noise=${noise} (already done) ==="
      continue
    fi
    echo "=== ca_only ${ckpt}, noise=${noise} ==="
    $PYBIN protein_mpnn_run.py --ca_only \
      --jsonl_path reproduction/phase1_data_prep/parsed_pdbs.jsonl \
      --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
      --out_folder "$outdir" \
      --model_name $ckpt \
      --num_seq_per_target 8 \
      --sampling_temp "0.1" \
      --seed 37 \
      --batch_size 1 \
      --backbone_noise $noise \
      --save_score 1 > "${outdir}.log" 2>&1
  done
done
echo "CHECKPOINT_SWEEP_ALL_DONE"

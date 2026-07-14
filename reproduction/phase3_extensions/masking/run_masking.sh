#!/bin/bash
set -e
cd /Users/christopherhuang/Documents/GitHub/ProteinMPNN
PYBIN=$(conda run -n proteinmpnn which python)
FRACS="10 20 30"

for frac in $FRACS; do
  echo "=== full-backbone, mask_frac=${frac}% ==="
  $PYBIN protein_mpnn_run.py \
    --jsonl_path reproduction/phase3_extensions/masking/masked_pdbs_frac${frac}.jsonl \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder reproduction/phase3_extensions/masking/fullbackbone_frac${frac} \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    --save_score 1 > reproduction/phase3_extensions/masking/fullbackbone_frac${frac}.log 2>&1

  echo "=== ca_only, mask_frac=${frac}% ==="
  $PYBIN protein_mpnn_run.py --ca_only \
    --jsonl_path reproduction/phase3_extensions/masking/masked_pdbs_frac${frac}.jsonl \
    --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
    --out_folder reproduction/phase3_extensions/masking/ca_only_frac${frac} \
    --model_name v_48_020 \
    --num_seq_per_target 8 \
    --sampling_temp "0.1" \
    --seed 37 \
    --batch_size 1 \
    --save_score 1 > reproduction/phase3_extensions/masking/ca_only_frac${frac}.log 2>&1
done

echo "=== baseline (0% mask, i.e. the original 21 clean proteins) for fair comparison ==="
conda run -n proteinmpnn python -c "
import json, csv
clean = set()
with open('reproduction/phase1_data_prep/dataset_summary.csv') as f:
    for r in csv.DictReader(f):
        if r['subset']=='clean': clean.add(r['pdb_id'])
with open('reproduction/phase1_data_prep/parsed_pdbs.jsonl') as f, open('reproduction/phase3_extensions/masking/masked_pdbs_frac0.jsonl','w') as out:
    for line in f:
        d = json.loads(line)
        if d['name'] in clean:
            out.write(line)
"
$PYBIN protein_mpnn_run.py \
  --jsonl_path reproduction/phase3_extensions/masking/masked_pdbs_frac0.jsonl \
  --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
  --out_folder reproduction/phase3_extensions/masking/fullbackbone_frac0 \
  --model_name v_48_020 --num_seq_per_target 8 --sampling_temp "0.1" --seed 37 --batch_size 1 --save_score 1 \
  > reproduction/phase3_extensions/masking/fullbackbone_frac0.log 2>&1
$PYBIN protein_mpnn_run.py --ca_only \
  --jsonl_path reproduction/phase3_extensions/masking/masked_pdbs_frac0.jsonl \
  --chain_id_jsonl reproduction/phase2_reproduction/chain_id.jsonl \
  --out_folder reproduction/phase3_extensions/masking/ca_only_frac0 \
  --model_name v_48_020 --num_seq_per_target 8 --sampling_temp "0.1" --seed 37 --batch_size 1 --save_score 1 \
  > reproduction/phase3_extensions/masking/ca_only_frac0.log 2>&1

echo "MASKING_ALL_DONE"

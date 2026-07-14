"""Parse Phase 2 FASTA outputs + terminal logs into per-sample metrics tables
for the full-backbone and CA-only models across the 0.1/0.2/0.3 temperature
sweep, plus a runtime table from the terminal logs."""
import re
import glob
import math
import csv
import os

MODELS = {
    "fullbackbone": "reproduction/phase2_reproduction/fullbackbone",
    "ca_only": "reproduction/phase2_reproduction/ca_only",
}

SAMPLE_RE = re.compile(
    r">T=([\d.]+), sample=(\d+), score=([\d.]+), global_score=([\d.]+), seq_recovery=([\d.]+)"
)
RUNTIME_RE = re.compile(r"^(\d+) sequences of length (\d+) generated in ([\d.]+) seconds")
PROTEIN_RE = re.compile(r"^Generating sequences for: (\S+)")

rows = []
for model_name, folder in MODELS.items():
    for fa_path in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
        pdb_id = os.path.basename(fa_path).replace(".fa", "")
        with open(fa_path) as f:
            for line in f:
                m = SAMPLE_RE.match(line.strip())
                if m:
                    temp, sample, score, global_score, seq_rec = m.groups()
                    score = float(score)
                    rows.append({
                        "model": model_name,
                        "pdb_id": pdb_id,
                        "temperature": float(temp),
                        "sample": int(sample),
                        "score_nll": score,
                        "global_score_nll": float(global_score),
                        "seq_recovery": float(seq_rec),
                        "perplexity": math.exp(score),
                    })

with open("reproduction/phase2_reproduction/per_sample_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {len(rows)} per-sample rows -> per_sample_metrics.csv")

# --- Runtime table from terminal logs ---
runtime_rows = []
for model_name, folder in MODELS.items():
    log_path = os.path.join(folder, "terminal_output.log")
    current_pdb = None
    with open(log_path) as f:
        for line in f:
            pm = PROTEIN_RE.match(line.strip())
            if pm:
                current_pdb = pm.group(1)
                continue
            rm = RUNTIME_RE.match(line.strip())
            if rm and current_pdb:
                n_seqs, length, seconds = rm.groups()
                n_seqs, length, seconds = int(n_seqs), int(length), float(seconds)
                runtime_rows.append({
                    "model": model_name,
                    "pdb_id": current_pdb,
                    "length": length,
                    "num_sequences": n_seqs,
                    "runtime_sec": seconds,
                    "seqs_per_sec": n_seqs / seconds if seconds > 0 else None,
                })
                current_pdb = None

with open("reproduction/phase2_reproduction/runtime_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(runtime_rows[0].keys()))
    w.writeheader()
    w.writerows(runtime_rows)
print(f"Wrote {len(runtime_rows)} runtime rows -> runtime_metrics.csv")

# --- Aggregate summary: mean seq_recovery / score / perplexity per model x temp ---
from collections import defaultdict
agg = defaultdict(list)
for r in rows:
    agg[(r["model"], r["temperature"])].append(r)

summary_rows = []
for (model, temp), items in sorted(agg.items()):
    n = len(items)
    summary_rows.append({
        "model": model,
        "temperature": temp,
        "n_samples": n,
        "mean_seq_recovery": sum(i["seq_recovery"] for i in items) / n,
        "mean_score_nll": sum(i["score_nll"] for i in items) / n,
        "mean_perplexity": sum(i["perplexity"] for i in items) / n,
    })

with open("reproduction/phase2_reproduction/temperature_sweep_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
    w.writeheader()
    w.writerows(summary_rows)
print("Wrote temperature_sweep_summary.csv")
for r in summary_rows:
    print(r)

# --- Runtime summary per model ---
rt_agg = defaultdict(list)
for r in runtime_rows:
    rt_agg[r["model"]].append(r)

print("\n--- Runtime summary ---")
for model, items in rt_agg.items():
    total_seqs = sum(i["num_sequences"] for i in items)
    total_time = sum(i["runtime_sec"] for i in items)
    print(f"{model}: total_proteins={len(items)} total_seqs={total_seqs} "
          f"total_time={total_time:.1f}s seqs_per_min={total_seqs/total_time*60:.1f} "
          f"mean_runtime_per_protein={total_time/len(items):.3f}s")

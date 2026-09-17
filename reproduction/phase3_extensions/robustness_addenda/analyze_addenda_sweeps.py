"""Analysis for the two follow-up sweeps that test whether the
noise-sweep full-backbone/CA-only crossover (~0.4-0.5 A) is real:

1. seed_sweep/   -- noise in {0.4, 0.45, 0.5, 0.6} x seed in {37, 38, 39},
   both models, at the matched v_48_020 checkpoint. Tests seed-stability of
   the crossover (the original sweep used a single seed=37 run per point).
2. checkpoint_sweep/ -- CA-only only, at v_48_002 and v_48_010 (different
   training noise levels than the matched v_48_020), noise in
   {0.0, 0.2, 0.3, 0.4, 0.5}. Tests whether the crossover point is an
   artifact of CA-only's specific v_48_020 training noise.

Reuses the same FASTA header parsing convention as
../noise_sweep/analyze_noise_sweep.py.
"""
import re
import glob
import csv
import os
from collections import defaultdict

SAMPLE_RE = re.compile(
    r">T=([\d.]+), sample=(\d+), score=([\d.]+), global_score=([\d.]+), seq_recovery=([\d.]+)"
)

HERE = os.path.dirname(os.path.abspath(__file__))
NOISE_SWEEP_SUMMARY = os.path.join(HERE, "..", "noise_sweep", "summary.csv")


def recoveries_for_folder(folder):
    vals = []
    for fa_path in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
        with open(fa_path) as f:
            for line in f:
                m = SAMPLE_RE.match(line.strip())
                if m:
                    vals.append(float(m.group(5)))
    return vals


# --- baseline (seed=37) recovery at each noise level, from the original sweep ---
baseline = {}  # (model, noise) -> mean recovery
with open(NOISE_SWEEP_SUMMARY) as f:
    for row in csv.DictReader(f):
        baseline[(row["model"], float(row["noise"]))] = float(row["mean_seq_recovery"])

# --- 1. seed sweep ---
SEED_NOISE_LEVELS = ["0.4", "0.45", "0.5", "0.6"]
SEEDS = [37, 38, 39]
seed_rows = []
for model in ["fullbackbone", "ca_only"]:
    for noise in SEED_NOISE_LEVELS:
        for seed in SEEDS:
            if seed == 37 and noise == "0.5":
                # reused from the original noise_sweep run (seed=37 default there)
                vals = None
                mean_rec = baseline[(model, 0.5)]
                n = None
            else:
                folder = os.path.join(HERE, "seed_sweep", f"{model}_noise{noise}_seed{seed}")
                vals = recoveries_for_folder(folder)
                mean_rec = sum(vals) / len(vals)
                n = len(vals)
            seed_rows.append({
                "model": model, "noise": float(noise), "seed": seed,
                "n_samples": n, "mean_seq_recovery": mean_rec,
            })

with open(os.path.join(HERE, "seed_sweep_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(seed_rows[0].keys()))
    w.writeheader()
    w.writerows(seed_rows)

# per-noise mean/std across seeds, and crossover check (fullbackbone - ca_only)
by_model_noise = defaultdict(list)
for r in seed_rows:
    by_model_noise[(r["model"], r["noise"])].append(r["mean_seq_recovery"])

crossover_lines = []
crossover_lines.append("| noise (A) | full-backbone mean (range across seeds) | CA-only mean (range across seeds) | FB - CA (pp) | which model wins |")
crossover_lines.append("|---|---|---|---|---|")
for noise in [float(n) for n in SEED_NOISE_LEVELS]:
    fb_vals = by_model_noise[("fullbackbone", noise)]
    ca_vals = by_model_noise[("ca_only", noise)]
    fb_mean = sum(fb_vals) / len(fb_vals)
    ca_mean = sum(ca_vals) / len(ca_vals)
    diff_pp = (fb_mean - ca_mean) * 100
    winner = "full-backbone" if diff_pp > 0 else "CA-only"
    crossover_lines.append(
        f"| {noise} | {fb_mean*100:.1f}% ({min(fb_vals)*100:.1f}-{max(fb_vals)*100:.1f}%) | "
        f"{ca_mean*100:.1f}% ({min(ca_vals)*100:.1f}-{max(ca_vals)*100:.1f}%) | {diff_pp:+.1f} | {winner} |"
    )

with open(os.path.join(HERE, "seed_sweep_crossover.md"), "w") as f:
    f.write("# Seed-stability of the noise-sweep crossover\n\n")
    f.write("Each cell is the mean over 3 seeds (37, 38, 39) x 30 structures x "
            "8 samples at T=0.1 (240 sequences per seed); noise=0.5/seed=37 "
            "reuses the original single-seed noise_sweep run rather than "
            "re-computing it.\n\n")
    f.write("\n".join(crossover_lines) + "\n")

print("=== Seed sweep crossover ===")
print("\n".join(crossover_lines))

# --- 2. checkpoint sweep (CA-only across training-noise checkpoints) ---
CKPT_NOISE_LEVELS = ["0.0", "0.2", "0.3", "0.4", "0.5"]
CHECKPOINTS = ["v_48_002", "v_48_010"]
ckpt_rows = []
for ckpt in CHECKPOINTS:
    for noise in CKPT_NOISE_LEVELS:
        folder = os.path.join(HERE, "checkpoint_sweep", f"ca_only_{ckpt}_noise{noise}")
        vals = recoveries_for_folder(folder)
        mean_rec = sum(vals) / len(vals)
        ckpt_rows.append({
            "checkpoint": ckpt, "noise": float(noise),
            "n_samples": len(vals), "mean_seq_recovery": mean_rec,
        })

with open(os.path.join(HERE, "checkpoint_sweep_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(ckpt_rows[0].keys()))
    w.writeheader()
    w.writerows(ckpt_rows)

ckpt_lines = []
ckpt_lines.append("| noise (A) | v_48_020 (matched, seed=37) | v_48_002 | v_48_010 |")
ckpt_lines.append("|---|---|---|---|")
for noise in [float(n) for n in CKPT_NOISE_LEVELS]:
    v020 = baseline.get(("ca_only", noise))
    v002 = next(r["mean_seq_recovery"] for r in ckpt_rows if r["checkpoint"] == "v_48_002" and r["noise"] == noise)
    v010 = next(r["mean_seq_recovery"] for r in ckpt_rows if r["checkpoint"] == "v_48_010" and r["noise"] == noise)
    v020_str = f"{v020*100:.1f}%" if v020 is not None else "n/a"
    ckpt_lines.append(f"| {noise} | {v020_str} | {v002*100:.1f}% | {v010*100:.1f}% |")

with open(os.path.join(HERE, "checkpoint_sweep_table.md"), "w") as f:
    f.write("# CA-only recovery across training-noise checkpoints\n\n")
    f.write("v_48_020/v_48_002/v_48_010 = CA-only checkpoints trained with "
            "Gaussian coordinate noise std. 0.20/0.02/0.10 A respectively. "
            "v_48_020 column reuses the original noise_sweep (single seed=37) "
            "run; v_48_002/v_48_010 use seed=37 only (this sweep's purpose is "
            "to see whether checkpoint choice shifts the curve, not to "
            "re-establish seed stability, which the seed_sweep already covers).\n\n")
    f.write("\n".join(ckpt_lines) + "\n")

print("\n=== Checkpoint sweep ===")
print("\n".join(ckpt_lines))
print("\nWrote seed_sweep_summary.csv, seed_sweep_crossover.md, "
      "checkpoint_sweep_summary.csv, checkpoint_sweep_table.md")

"""Phase 3 robustness: parse noise-sweep FASTA outputs into a metrics table
and plot recovery/perplexity degradation vs. backbone noise level, comparing
full-backbone and CA-only models."""
import re
import glob
import math
import csv
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NOISE_LEVELS = ["0.0", "0.1", "0.2", "0.3", "0.5"]
SAMPLE_RE = re.compile(
    r">T=([\d.]+), sample=(\d+), score=([\d.]+), global_score=([\d.]+), seq_recovery=([\d.]+)"
)

rows = []
for model in ["fullbackbone", "ca_only"]:
    for noise in NOISE_LEVELS:
        folder = f"reproduction/phase3_extensions/noise_sweep/{model}_noise{noise}"
        for fa_path in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
            pdb_id = os.path.basename(fa_path).replace(".fa", "")
            with open(fa_path) as f:
                for line in f:
                    m = SAMPLE_RE.match(line.strip())
                    if m:
                        _, sample, score, global_score, seq_rec = m.groups()
                        score = float(score)
                        rows.append({
                            "model": model, "noise": float(noise), "pdb_id": pdb_id,
                            "sample": int(sample), "score_nll": score,
                            "seq_recovery": float(seq_rec), "perplexity": math.exp(score),
                        })

with open("reproduction/phase3_extensions/noise_sweep/per_sample_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {len(rows)} rows -> per_sample_metrics.csv")

agg = defaultdict(list)
for r in rows:
    agg[(r["model"], r["noise"])].append(r)

summary = []
for (model, noise), items in sorted(agg.items(), key=lambda x: (x[0][0], x[0][1])):
    n = len(items)
    summary.append({
        "model": model, "noise": noise, "n_samples": n,
        "mean_seq_recovery": sum(i["seq_recovery"] for i in items) / n,
        "mean_perplexity": sum(i["perplexity"] for i in items) / n,
    })
with open("reproduction/phase3_extensions/noise_sweep/summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader()
    w.writerows(summary)
for r in summary:
    print(r)

COLOR = {"fullbackbone": "#2a78d6", "ca_only": "#1baf7a"}
LABEL = {"fullbackbone": "Full-backbone", "ca_only": "CA-only"}
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
for model in ["fullbackbone", "ca_only"]:
    items = [s for s in summary if s["model"] == model]
    xs = [s["noise"] for s in items]
    axes[0].plot(xs, [s["mean_seq_recovery"] * 100 for s in items], marker="o",
                 markersize=7, linewidth=2, color=COLOR[model], label=LABEL[model])
    axes[1].plot(xs, [s["mean_perplexity"] for s in items], marker="o",
                 markersize=7, linewidth=2, color=COLOR[model], label=LABEL[model])
axes[0].set_xlabel("Gaussian backbone noise (Å std. dev.)")
axes[0].set_ylabel("Mean sequence recovery (%)")
axes[0].set_title("Recovery vs. backbone noise")
axes[1].set_xlabel("Gaussian backbone noise (Å std. dev.)")
axes[1].set_ylabel("Mean perplexity")
axes[1].set_title("Perplexity vs. backbone noise")
for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
    ax.set_axisbelow(True)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2, frameon=False)
fig.tight_layout()
fig.savefig("reproduction/phase3_extensions/noise_sweep/noise_sweep_comparison.png", dpi=200, bbox_inches="tight")
print("Wrote noise_sweep_comparison.png")

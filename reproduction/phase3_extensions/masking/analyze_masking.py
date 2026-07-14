"""Phase 3 robustness: parse masking-experiment FASTA outputs, plot recovery
degradation vs. mask fraction for full-backbone vs. CA-only. Recovery is
computed by the official pipeline only over non-masked (isfinite-coordinate)
positions -- see protein_mpnn_utils.py:_scores and tied_featurize's `mask`
variable -- so this directly measures "how well does the model design the
still-visible part of the structure as more of the rest is hidden"."""
import re
import glob
import math
import csv
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FRACS = ["0", "10", "20", "30"]
SAMPLE_RE = re.compile(
    r">T=([\d.]+), sample=(\d+), score=([\d.]+), global_score=([\d.]+), seq_recovery=([\d.]+)"
)

rows = []
for model in ["fullbackbone", "ca_only"]:
    for frac in FRACS:
        folder = f"reproduction/phase3_extensions/masking/{model}_frac{frac}"
        for fa_path in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
            pdb_id = os.path.basename(fa_path).replace(".fa", "")
            with open(fa_path) as f:
                for line in f:
                    m = SAMPLE_RE.match(line.strip())
                    if m:
                        _, sample, score, global_score, seq_rec = m.groups()
                        score = float(score)
                        rows.append({
                            "model": model, "mask_frac": int(frac), "pdb_id": pdb_id,
                            "sample": int(sample), "score_nll": score,
                            "seq_recovery": float(seq_rec), "perplexity": math.exp(score),
                        })

with open("reproduction/phase3_extensions/masking/per_sample_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {len(rows)} rows -> per_sample_metrics.csv")

agg = defaultdict(list)
for r in rows:
    agg[(r["model"], r["mask_frac"])].append(r)

summary = []
for (model, frac), items in sorted(agg.items(), key=lambda x: (x[0][0], x[0][1])):
    n = len(items)
    summary.append({
        "model": model, "mask_frac_pct": frac, "n_samples": n,
        "n_proteins": n // 8,
        "mean_seq_recovery": sum(i["seq_recovery"] for i in items) / n,
        "mean_perplexity": sum(i["perplexity"] for i in items) / n,
    })
with open("reproduction/phase3_extensions/masking/summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader()
    w.writerows(summary)
for r in summary:
    print(r)

COLOR = {"fullbackbone": "#2a78d6", "ca_only": "#1baf7a"}
LABEL = {"fullbackbone": "Full-backbone", "ca_only": "CA-only"}
fig, ax = plt.subplots(figsize=(7, 4.2))
for model in ["fullbackbone", "ca_only"]:
    items = [s for s in summary if s["model"] == model]
    xs = [s["mask_frac_pct"] for s in items]
    ax.plot(xs, [s["mean_seq_recovery"] * 100 for s in items], marker="o",
            markersize=7, linewidth=2, color=COLOR[model], label=LABEL[model])
ax.set_xlabel("Residues masked (%)")
ax.set_ylabel("Mean recovery on remaining visible positions (%)")
ax.set_title("Recovery vs. synthetic residue masking (21 clean proteins)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
ax.set_axisbelow(True)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("reproduction/phase3_extensions/masking/masking_comparison.png", dpi=200, bbox_inches="tight")
print("Wrote masking_comparison.png")

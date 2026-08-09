"""Phase 2 comparison target (2nd half): does model behavior change
predictably with input completeness? Joins per_sample_metrics.csv against
the Phase 1 clean/stress_test subset labels."""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

subset = {}
with open("reproduction/phase1_data_prep/dataset_summary.csv") as f:
    for r in csv.DictReader(f):
        subset[r["pdb_id"]] = r["subset"]

rows = list(csv.DictReader(open("reproduction/phase2_reproduction/per_sample_metrics.csv")))
for r in rows:
    r["subset"] = subset[r["pdb_id"]]

agg = defaultdict(list)
for r in rows:
    agg[(r["model"], r["temperature"], r["subset"])].append(r)

summary_rows = []
for (model, temp, sub), items in sorted(agg.items(), key=lambda x: (x[0][0], float(x[0][1]), x[0][2])):
    n = len(items)
    summary_rows.append({
        "model": model,
        "temperature": temp,
        "subset": sub,
        "n_samples": n,
        "n_proteins": n // 8,
        "mean_seq_recovery": sum(float(i["seq_recovery"]) for i in items) / n,
        "mean_score_nll": sum(float(i["score_nll"]) for i in items) / n,
        "mean_perplexity": sum(float(i["perplexity"]) for i in items) / n,
    })

with open("reproduction/phase2_reproduction/completeness_comparison.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
    w.writeheader()
    w.writerows(summary_rows)
print("Wrote completeness_comparison.csv")
for r in summary_rows:
    print(r)

# --- Plot: mean recovery, clean vs stress_test, per model, at each temperature ---
COLOR = {"clean": "#2a78d6", "stress_test": "#e34948"}
LABEL = {"clean": "Clean (21 proteins)", "stress_test": "Stress-test (9 proteins)"}
MODEL_LABEL = {"fullbackbone": "Full-backbone", "ca_only": "CA-only"}

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
temps = [0.1, 0.2, 0.3]
width = 0.35
for ax, model in zip(axes, ["fullbackbone", "ca_only"]):
    for i, sub in enumerate(["clean", "stress_test"]):
        ys = [next(r["mean_seq_recovery"] for r in summary_rows
                   if r["model"] == model and float(r["temperature"]) == t and r["subset"] == sub)
              for t in temps]
        xs = [j + (i - 0.5) * width for j in range(len(temps))]
        ax.bar(xs, [y * 100 for y in ys], width=width, color=COLOR[sub], label=LABEL[sub], zorder=3)
    ax.set_xticks(range(len(temps)))
    ax.set_xticklabels([str(t) for t in temps])
    ax.set_xlabel("Sampling temperature")
    ax.set_title(MODEL_LABEL[model])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
    ax.set_axisbelow(True)
axes[0].set_ylabel("Mean sequence recovery (%)")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=2, frameon=False)
fig.suptitle("")
fig.tight_layout()
fig.savefig("reproduction/phase2_reproduction/completeness_comparison.png", dpi=200, bbox_inches="tight")
print("Wrote completeness_comparison.png")

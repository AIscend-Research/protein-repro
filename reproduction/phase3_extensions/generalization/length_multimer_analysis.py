"""Phase 3 generalization analysis: recovery vs. protein length and
monomer vs. multimer, reusing Phase 2's per_sample_metrics.csv (no new
inference runs needed) joined against Phase 1's dataset_summary.csv."""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

meta = {}
with open("reproduction/phase1_data_prep/dataset_summary.csv") as f:
    for r in csv.DictReader(f):
        meta[r["pdb_id"]] = r

rows = list(csv.DictReader(open("reproduction/phase2_reproduction/per_sample_metrics.csv")))
# Use T=0.1 (standard temperature) only, to isolate the length/multimer effect
rows = [r for r in rows if r["temperature"] == "0.1"]

for r in rows:
    m = meta[r["pdb_id"]]
    r["length"] = int(m["total_resolved_residues"])
    r["monomer_or_multimer"] = m["monomer_or_multimer"]

LENGTH_BINS = [(0, 70, "short (<70)"), (70, 150, "medium (70-150)"), (150, 10**6, "long (>150)")]


def length_bucket(L):
    for lo, hi, label in LENGTH_BINS:
        if lo <= L < hi:
            return label
    return "?"


for r in rows:
    r["length_bucket"] = length_bucket(r["length"])

print("=== Recovery by length bucket (T=0.1) ===")
by_bucket = defaultdict(list)
for r in rows:
    by_bucket[(r["model"], r["length_bucket"])].append(float(r["seq_recovery"]))

bucket_order = [b[2] for b in LENGTH_BINS]
length_summary = []
for model in ["fullbackbone", "ca_only"]:
    for bucket in bucket_order:
        vals = by_bucket[(model, bucket)]
        n_proteins = len(vals) // 8
        mean_r = sum(vals) / len(vals)
        length_summary.append({"model": model, "length_bucket": bucket,
                                "n_proteins": n_proteins, "mean_seq_recovery": mean_r})
        print(f"{model:14}{bucket:20}n_proteins={n_proteins:<4}mean_recovery={mean_r:.4f}")

with open("reproduction/phase3_extensions/generalization/length_bucket_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(length_summary[0].keys()))
    w.writeheader()
    w.writerows(length_summary)

print("\n=== Recovery by monomer vs multimer (T=0.1) ===")
by_mm = defaultdict(list)
for r in rows:
    by_mm[(r["model"], r["monomer_or_multimer"])].append(float(r["seq_recovery"]))

mm_summary = []
for model in ["fullbackbone", "ca_only"]:
    for mm in ["monomer", "multimer"]:
        vals = by_mm[(model, mm)]
        n_proteins = len(vals) // 8
        mean_r = sum(vals) / len(vals)
        mm_summary.append({"model": model, "monomer_or_multimer": mm,
                            "n_proteins": n_proteins, "mean_seq_recovery": mean_r})
        print(f"{model:14}{mm:14}n_proteins={n_proteins:<4}mean_recovery={mean_r:.4f}")

with open("reproduction/phase3_extensions/generalization/monomer_multimer_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(mm_summary[0].keys()))
    w.writeheader()
    w.writerows(mm_summary)

# --- Plot: length bucket comparison ---
COLOR = {"fullbackbone": "#2a78d6", "ca_only": "#1baf7a"}
fig, ax = plt.subplots(figsize=(7, 4.2))
x = list(range(len(bucket_order)))
width = 0.35
for i, model in enumerate(["fullbackbone", "ca_only"]):
    ys = [next(s["mean_seq_recovery"] for s in length_summary
               if s["model"] == model and s["length_bucket"] == b) * 100 for b in bucket_order]
    xs = [j + (i - 0.5) * width for j in x]
    ax.bar(xs, ys, width=width, color=COLOR[model],
           label="Full-backbone" if model == "fullbackbone" else "CA-only", zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(bucket_order)
ax.set_ylabel("Mean sequence recovery (%)")
ax.set_xlabel("Protein length bucket")
ax.set_title("Recovery vs. protein length (T=0.1)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
ax.set_axisbelow(True)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("reproduction/phase3_extensions/generalization/length_bucket_comparison.png", dpi=200, bbox_inches="tight")
print("\nWrote length_bucket_comparison.png")

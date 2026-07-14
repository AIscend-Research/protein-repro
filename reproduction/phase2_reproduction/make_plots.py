"""Phase 2 comparison plots: full-backbone vs CA-only, sequence recovery and
perplexity across the temperature sweep. Two single-axis subplots (never a
dual-axis chart), fixed categorical colors (blue=full-backbone, aqua=CA-only)."""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLOR = {"fullbackbone": "#2a78d6", "ca_only": "#1baf7a"}
LABEL = {"fullbackbone": "Full-backbone (v_48_020)", "ca_only": "CA-only (v_48_020)"}

rows = list(csv.DictReader(open("reproduction/phase2_reproduction/temperature_sweep_summary.csv")))
by_model = defaultdict(list)
for r in rows:
    by_model[r["model"]].append(r)
for m in by_model:
    by_model[m].sort(key=lambda r: float(r["temperature"]))

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

ax = axes[0]
for model in ["fullbackbone", "ca_only"]:
    xs = [float(r["temperature"]) for r in by_model[model]]
    ys = [float(r["mean_seq_recovery"]) * 100 for r in by_model[model]]
    ax.plot(xs, ys, marker="o", markersize=7, linewidth=2, color=COLOR[model], label=LABEL[model])
ax.set_xlabel("Sampling temperature")
ax.set_ylabel("Mean sequence recovery (%)")
ax.set_title("Sequence recovery vs. temperature")
ax.set_xticks([0.1, 0.2, 0.3])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
ax.set_axisbelow(True)

ax = axes[1]
for model in ["fullbackbone", "ca_only"]:
    xs = [float(r["temperature"]) for r in by_model[model]]
    ys = [float(r["mean_perplexity"]) for r in by_model[model]]
    ax.plot(xs, ys, marker="o", markersize=7, linewidth=2, color=COLOR[model], label=LABEL[model])
ax.set_xlabel("Sampling temperature")
ax.set_ylabel("Mean perplexity")
ax.set_title("Perplexity vs. temperature")
ax.set_xticks([0.1, 0.2, 0.3])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
ax.set_axisbelow(True)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=2, frameon=False)
fig.tight_layout()
fig.savefig("reproduction/phase2_reproduction/temperature_sweep_comparison.png", dpi=200, bbox_inches="tight")
print("Wrote temperature_sweep_comparison.png")

# --- Runtime bar comparison ---
rt_rows = list(csv.DictReader(open("reproduction/phase2_reproduction/runtime_metrics.csv")))
rt_by_model = defaultdict(list)
for r in rt_rows:
    rt_by_model[r["model"]].append(r)

fig2, ax2 = plt.subplots(figsize=(6, 4.2))
models = ["fullbackbone", "ca_only"]
means = [sum(float(r["seqs_per_sec"]) for r in rt_by_model[m]) / len(rt_by_model[m]) for m in models]
bars = ax2.bar([LABEL[m] for m in models], means, color=[COLOR[m] for m in models], width=0.5, zorder=3)
ax2.set_ylabel("Mean sequences / second / protein")
ax2.set_title("Inference throughput (CPU/MPS, no CUDA)")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
ax2.set_axisbelow(True)
for b, v in zip(bars, means):
    ax2.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}", ha="center", va="bottom")
fig2.tight_layout()
fig2.savefig("reproduction/phase2_reproduction/runtime_comparison.png", dpi=200, bbox_inches="tight")
print("Wrote runtime_comparison.png")

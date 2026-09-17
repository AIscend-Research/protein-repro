"""Phase 3 low-resource inference: parse the run_low_resource.sh logs into
runtime and peak-memory tables, and plot the two low-compute levers
(num_seq_per_target and batch_size) against wall time.

Reads the `.log` files written by run_low_resource.sh (plus the Phase 2
memory_check logs, if present) and writes:
  runtime_summary.csv  - per-condition total/mean inference time
  memory_summary.csv   - peak RSS / memory footprint per condition
  low_resource_comparison.png
"""
import csv
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "reproduction/phase3_extensions/low_resource"
MEM_BASE = "reproduction/phase2_reproduction/memory_check"

# "8 sequences of length 76 generated in 1.1423 seconds"
GEN_RE = re.compile(r"^(\d+) sequences of length (\d+) generated in ([\d.]+) seconds")
# bash `time` builtin: "real\t0m5.720s"
REAL_RE = re.compile(r"^real\s+(\d+)m([\d.]+)s")
# /usr/bin/time -l (macOS)
RSS_RE = re.compile(r"^\s*(\d+)\s+maximum resident set size")
FOOTPRINT_RE = re.compile(r"^\s*(\d+)\s+peak memory footprint")

# label -> (condition group, num_seq_per_target, batch_size)
RUNTIME_CASES = {
    "highcompute_nseq8": ("num_seq", 8, 1),
    "lowcompute_nseq1": ("num_seq", 1, 1),
    "batchsize1": ("batch_size", 8, 1),
    "batchsize2": ("batch_size", 8, 2),
    "batchsize4": ("batch_size", 8, 4),
    "batchsize8": ("batch_size", 8, 8),
}

MEMORY_CASES = {
    # label -> (log path, structure, num_seq, batch_size, model)
    "mem_bs1": (f"{BASE}/mem_bs1.log", "3PGK", 8, 1, "fullbackbone"),
    "mem_bs8": (f"{BASE}/mem_bs8.log", "3PGK", 8, 8, "fullbackbone"),
    "fullbackbone_3PGK": (f"{MEM_BASE}/fullbackbone_3PGK.log", "3PGK", 8, 1, "fullbackbone"),
    "ca_only_3PGK": (f"{MEM_BASE}/ca_only_3PGK.log", "3PGK", 8, 1, "ca_only"),
    "fullbackbone_1TUP": (f"{MEM_BASE}/fullbackbone_1TUP.log", "1TUP", 8, 1, "fullbackbone"),
}


def parse_runtime(path):
    """-> (summed per-protein inference seconds, n proteins, wall-clock seconds)."""
    per_protein, wall = [], None
    with open(path) as f:
        for line in f:
            m = GEN_RE.match(line.strip())
            if m:
                per_protein.append(float(m.group(3)))
                continue
            m = REAL_RE.match(line.strip())
            if m:
                wall = int(m.group(1)) * 60 + float(m.group(2))
    return sum(per_protein), len(per_protein), wall


def parse_memory(path):
    """-> (max RSS bytes, peak footprint bytes); either may be None."""
    rss = footprint = None
    with open(path) as f:
        for line in f:
            m = RSS_RE.match(line)
            if m:
                rss = int(m.group(1))
            m = FOOTPRINT_RE.match(line)
            if m:
                footprint = int(m.group(1))
    return rss, footprint


missing = []

runtime_rows = []
for label, (group, num_seq, batch) in RUNTIME_CASES.items():
    path = f"{BASE}/{label}.log"
    if not os.path.exists(path):
        missing.append(path)
        continue
    total, n, wall = parse_runtime(path)
    runtime_rows.append({
        "label": label,
        "lever": group,
        "num_seq_per_target": num_seq,
        "batch_size": batch,
        "n_proteins": n,
        "total_inference_s": round(total, 3),
        "mean_s_per_protein": round(total / n, 3) if n else None,
        "wall_clock_s": round(wall, 3) if wall is not None else None,
        "seqs_per_min": round(num_seq * n / total * 60, 1) if total else None,
    })

memory_rows = []
for label, (path, structure, num_seq, batch, model) in MEMORY_CASES.items():
    if not os.path.exists(path):
        missing.append(path)
        continue
    rss, footprint = parse_memory(path)
    memory_rows.append({
        "label": label,
        "model": model,
        "structure": structure,
        "num_seq_per_target": num_seq,
        "batch_size": batch,
        "max_rss_mb": round(rss / 1e6, 1) if rss else None,
        "peak_footprint_mb": round(footprint / 1e6, 1) if footprint else None,
    })


def write_csv(path, rows):
    if not rows:
        print(f"!! no rows for {path} — skipped")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows -> {path}")
    for r in rows:
        print("  ", r)


write_csv(f"{BASE}/runtime_summary.csv", runtime_rows)
write_csv(f"{BASE}/memory_summary.csv", memory_rows)

if missing:
    print("\n!! missing logs (run run_low_resource.sh / run_memory_check.sh first):")
    for p in missing:
        print("   ", p)

# --- plot: the two levers side by side ---------------------------------------
batch_rows = sorted((r for r in runtime_rows if r["lever"] == "batch_size"),
                    key=lambda r: r["batch_size"])
nseq_rows = sorted((r for r in runtime_rows if r["lever"] == "num_seq"),
                   key=lambda r: r["num_seq_per_target"])
mem_by_label = {r["label"]: r for r in memory_rows}

if batch_rows:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

    xs = [r["batch_size"] for r in batch_rows]
    axes[0].plot(xs, [r["total_inference_s"] for r in batch_rows], marker="o",
                 markersize=7, linewidth=2, color="#2a78d6", label="Inference time (s)")
    axes[0].set_xlabel("--batch_size (num_seq_per_target = 8)")
    axes[0].set_ylabel("Total time, 5 proteins (s)")
    axes[0].set_title("Speed vs. batch size")
    axes[0].set_xticks(xs)

    # Memory on the same panel, second axis — the tradeoff is the whole point.
    mem_xs, mem_ys = [], []
    for bs in (1, 8):
        row = mem_by_label.get(f"mem_bs{bs}")
        if row and row["max_rss_mb"]:
            mem_xs.append(bs)
            mem_ys.append(row["max_rss_mb"])
    if mem_ys:
        ax2 = axes[0].twinx()
        ax2.plot(mem_xs, mem_ys, marker="s", markersize=7, linewidth=2,
                 linestyle="--", color="#d6602a", label="Peak RSS (MB)")
        ax2.set_ylabel("Peak RSS on 3PGK (MB)")
        ax2.legend(loc="upper center", frameon=False)
    axes[0].legend(loc="upper right", frameon=False)

    if nseq_rows:
        labels = [str(r["num_seq_per_target"]) for r in nseq_rows]
        axes[1].bar(labels, [r["total_inference_s"] for r in nseq_rows],
                    color="#1baf7a", width=0.55)
        for i, r in enumerate(nseq_rows):
            axes[1].text(i, r["total_inference_s"], f' {r["total_inference_s"]:.2f}s',
                         ha="center", va="bottom")
    axes[1].set_xlabel("--num_seq_per_target (batch_size = 1)")
    axes[1].set_ylabel("Total time, 5 proteins (s)")
    axes[1].set_title("Speed vs. sequences per backbone")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.grid(axis="y", color="#e3e2dd", linewidth=1, zorder=0)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{BASE}/low_resource_comparison.png", dpi=200, bbox_inches="tight")
    print("Wrote low_resource_comparison.png")
else:
    print("!! no batch-size runtime rows — plot skipped")

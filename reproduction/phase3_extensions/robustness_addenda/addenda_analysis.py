"""Reviewer-driven addenda (no new inference needed for this script; all
numbers come from existing Phase 2 / Phase 3 per-sample CSVs and noise-sweep
FASTA outputs):

1. Paired per-structure recovery differences (full-backbone vs. CA-only),
   with win counts and a sign test, for the T=0.1 temperature-sweep headline
   gap and for each noise level of the noise sweep. Replaces "mean A - mean
   B" with a distribution across the 30 structures.
2. Duplicate-identity-corrected headline recovery: 1HEL (identical to 1AKI)
   and 1UBQ (identical to 1UBI) dropped, so the "30-structure" mean becomes
   a genuine 28-unique-sequence mean.
3. Length-reweighted recovery gap bound: recovery restricted to the
   medium-length (70-150 aa) bucket, the closest proxy available to a
   typical CATH-domain-sized target, to bound how much of the 6-point gap
   to the published 52.4% is attributable to length-distribution mismatch
   rather than genuine dataset mismatch.
4. Amino-acid composition Shannon entropy per (model, noise level), to check
   whether the 0.4-0.5 A noise crossover reflects a real robustness
   difference or both models collapsing toward degenerate low-entropy
   compositions.
"""
import csv
import glob
import math
import os
import re
from collections import Counter, defaultdict

REPO = "/Users/christopherhuang/Documents/GitHub/mlrc-ProteinMPNN"
OUT = os.path.join(REPO, "reproduction/phase3_extensions/robustness_addenda")
DUPLICATES = {"1HEL": "1AKI", "1UBQ": "1UBI"}  # drop key, identical to value

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")


def sign_test_p(n_pos, n_neg):
    """Two-sided exact binomial sign test p-value, p=0.5, no scipy dependency."""
    n = n_pos + n_neg
    if n == 0:
        return 1.0
    k = min(n_pos, n_neg)

    def choose(n, k):
        return math.comb(n, k)

    tail = sum(choose(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def paired_diff_report(rows, group_key, out_path, header_note):
    """rows: list of dicts with model, pdb_id, seq_recovery, plus group_key field."""
    by_group_model_pdb = defaultdict(dict)
    for r in rows:
        g = r[group_key]
        by_group_model_pdb[(g, r["model"])].setdefault(r["pdb_id"], []).append(r["seq_recovery"])

    groups = sorted(set(r[group_key] for r in rows))
    lines = [header_note, ""]
    lines.append(f"| {group_key} | n_pdb | mean(FB-CA) pp | median(FB-CA) pp | FB wins | CA wins | ties | sign-test p |")
    lines.append("|---|---|---|---|---|---|---|---|")
    summary_rows = []
    for g in groups:
        fb = by_group_model_pdb.get((g, "fullbackbone")) or by_group_model_pdb.get((g, "full_backbone"))
        ca = by_group_model_pdb.get((g, "ca_only"))
        if not fb or not ca:
            continue
        pdbs = sorted(set(fb) & set(ca))
        diffs = []
        for p in pdbs:
            fb_mean = sum(fb[p]) / len(fb[p])
            ca_mean = sum(ca[p]) / len(ca[p])
            diffs.append(fb_mean - ca_mean)
        n_pos = sum(1 for d in diffs if d > 0)
        n_neg = sum(1 for d in diffs if d < 0)
        n_tie = sum(1 for d in diffs if d == 0)
        p_val = sign_test_p(n_pos, n_neg)
        mean_d = sum(diffs) / len(diffs) * 100
        sorted_d = sorted(diffs)
        med_d = sorted_d[len(sorted_d) // 2] * 100
        lines.append(f"| {g} | {len(pdbs)} | {mean_d:+.1f} | {med_d:+.1f} | {n_pos} | {n_neg} | {n_tie} | {p_val:.4f} |")
        summary_rows.append({"group": g, "n_pdb": len(pdbs), "mean_diff_pp": mean_d,
                              "median_diff_pp": med_d, "fb_wins": n_pos, "ca_wins": n_neg,
                              "ties": n_tie, "sign_test_p": p_val})
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return summary_rows


# ---------------------------------------------------------------------------
# 1+2. Phase 2 temperature sweep: paired diffs + dedup headline recovery
# ---------------------------------------------------------------------------
phase2_rows = []
with open(os.path.join(REPO, "reproduction/phase2_reproduction/per_sample_metrics.csv")) as f:
    for r in csv.DictReader(f):
        r["temperature"] = r["temperature"]
        r["seq_recovery"] = float(r["seq_recovery"])
        phase2_rows.append(r)

t01_rows = [r for r in phase2_rows if r["temperature"] == "0.1"]
paired_diff_report(
    t01_rows, "temperature",
    os.path.join(OUT, "phase2_paired_diffs.md"),
    "# Phase 2 temperature sweep: paired per-structure full-backbone vs. CA-only recovery\n"
    "(FB-CA = full-backbone minus CA-only mean recovery per structure, percentage points)",
)
# same, but across all 3 temperatures for completeness
all_temp_summary = paired_diff_report(
    phase2_rows, "temperature",
    os.path.join(OUT, "phase2_paired_diffs_all_temps.md"),
    "# Phase 2, all temperatures: paired per-structure full-backbone vs. CA-only recovery",
)

# dedup headline recovery at T=0.1
dedup_rows = [r for r in t01_rows if r["pdb_id"] not in DUPLICATES]
headline = {}
for label, rowset in [("all_30", t01_rows), ("dedup_28", dedup_rows)]:
    agg = defaultdict(list)
    for r in rowset:
        agg[r["model"]].append(r["seq_recovery"])
    headline[label] = {m: sum(v) / len(v) * 100 for m, v in agg.items()}

with open(os.path.join(OUT, "dedup_headline_recovery.md"), "w") as f:
    f.write("# Duplicate-identity-corrected headline recovery (T=0.1)\n\n")
    f.write("1HEL is sequence-identical to 1AKI; 1UBQ is sequence-identical to 1UBI. "
            "\"dedup_28\" drops 1HEL and 1UBQ so each unique protein is counted once.\n\n")
    f.write("| subset | full-backbone | ca_only |\n|---|---|---|\n")
    for label in ["all_30", "dedup_28"]:
        fb = headline[label].get("fullbackbone", float("nan"))
        ca = headline[label].get("ca_only", float("nan"))
        f.write(f"| {label} | {fb:.1f}% | {ca:.1f}% |\n")

# ---------------------------------------------------------------------------
# 3. Length-reweighted recovery gap bound (medium-length bucket only)
# ---------------------------------------------------------------------------
lengths = {}
with open(os.path.join(REPO, "reproduction/phase1_data_prep/dataset_summary.csv")) as f:
    for r in csv.DictReader(f):
        lengths[r["pdb_id"]] = int(r["total_resolved_residues"])


def bucket(n):
    if n < 70:
        return "short"
    if n <= 150:
        return "medium"
    return "long"


bucketed = defaultdict(list)
for r in t01_rows:
    L = lengths.get(r["pdb_id"])
    if L is None:
        continue
    bucketed[(r["model"], bucket(L))].append(r["seq_recovery"])

PUBLISHED_RECOVERY = 52.4
with open(os.path.join(OUT, "length_reweighted_gap.md"), "w") as f:
    f.write("# Length-reweighted recovery gap bound\n\n")
    f.write(f"Published CATH-test-split recovery: {PUBLISHED_RECOVERY}%. "
            "Medium-length (70-150 aa) structures are the closest available proxy to "
            "typical CATH-domain sizes; comparing recovery on that subset alone to the "
            "full 30-structure mean bounds how much of the gap is length-distribution "
            "skew rather than genuine train/test mismatch.\n\n")
    f.write("| model | all-30 mean | medium-only mean | shift toward published | gap all-30 | gap medium-only |\n")
    f.write("|---|---|---|---|---|---|\n")
    for model in ["fullbackbone", "ca_only"]:
        all_vals = [r["seq_recovery"] for r in t01_rows if r["model"] == model]
        med_vals = bucketed[(model, "medium")]
        all_mean = sum(all_vals) / len(all_vals) * 100
        med_mean = sum(med_vals) / len(med_vals) * 100
        shift = med_mean - all_mean
        gap_all = PUBLISHED_RECOVERY - all_mean
        gap_med = PUBLISHED_RECOVERY - med_mean
        f.write(f"| {model} | {all_mean:.1f}% | {med_mean:.1f}% | {shift:+.1f} pp | "
                f"{gap_all:+.1f} pp | {gap_med:+.1f} pp |\n")
    f.write(f"\n(full-backbone is the model paper's headline comparison; "
            f"published number is full-backbone T=0.1.)\n")

# ---------------------------------------------------------------------------
# 4. AA composition entropy per (model, noise) from noise-sweep FASTAs
# ---------------------------------------------------------------------------
SAMPLE_RE = re.compile(r">T=([\d.]+), sample=(\d+),")
NOISE_LEVELS = ["0.0", "0.1", "0.2", "0.3", "0.5"]


def shannon_entropy(seq):
    counts = Counter(seq)
    n = len(seq)
    ent = 0.0
    for c in counts.values():
        p = c / n
        ent -= p * math.log2(p)
    return ent


MAX_ENTROPY = math.log2(20)
entropy_rows = []
for model in ["fullbackbone", "ca_only"]:
    for noise in NOISE_LEVELS:
        folder = os.path.join(REPO, f"reproduction/phase3_extensions/noise_sweep/{model}_noise{noise}/seqs")
        for fa_path in sorted(glob.glob(os.path.join(folder, "*.fa"))):
            pdb_id = os.path.basename(fa_path).replace(".fa", "")
            with open(fa_path) as f:
                lines = f.readlines()
            for i in range(len(lines)):
                if lines[i].startswith(">T="):
                    seq = lines[i + 1].strip()
                    entropy_rows.append({
                        "model": model, "noise": float(noise), "pdb_id": pdb_id,
                        "entropy_bits": shannon_entropy(seq),
                        "normalized_entropy": shannon_entropy(seq) / MAX_ENTROPY,
                    })

agg_entropy = defaultdict(list)
for r in entropy_rows:
    agg_entropy[(r["model"], r["noise"])].append(r["normalized_entropy"])

with open(os.path.join(OUT, "noise_sweep_entropy.md"), "w") as f:
    f.write("# Amino-acid composition entropy vs. backbone noise\n\n")
    f.write("Normalized Shannon entropy (0-1, 1 = uniform over 20 AAs) of designed-sequence "
            "composition, averaged per (model, noise level). Tests whether the 0.4-0.5 A "
            "recovery crossover coincides with a composition collapse.\n\n")
    f.write("| model | noise (A) | mean normalized entropy |\n|---|---|---|\n")
    for (model, noise), vals in sorted(agg_entropy.items(), key=lambda x: (x[0][0], x[0][1])):
        f.write(f"| {model} | {noise} | {sum(vals)/len(vals):.4f} |\n")

with open(os.path.join(OUT, "noise_sweep_entropy_per_sample.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["model", "noise", "pdb_id", "entropy_bits", "normalized_entropy"])
    w.writeheader()
    w.writerows(entropy_rows)

print("Wrote paired-diff, dedup, length-reweight, and entropy addenda to", OUT)
for fn in sorted(os.listdir(OUT)):
    print(" -", fn)

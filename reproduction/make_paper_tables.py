"""Regenerate every table in the write-up from the result CSVs, so no number in
the paper is transcribed by hand.

Reads the summary CSVs produced by the per-phase analysis scripts and writes
`reproduction/paper_tables.md` — one markdown table per paper table, each
captioned with its source file so a reader can trace any figure back to the
run that produced it. Missing inputs are reported and skipped rather than
faked, so a partial run still produces a usable (and honest) document.

Usage: python reproduction/make_paper_tables.py
"""
import csv
import os
from collections import OrderedDict

OUT_PATH = "reproduction/paper_tables.md"

MODEL_LABEL = {"fullbackbone": "Full-backbone", "ca_only": "CA-only"}
# Columns rendered as percentages rather than raw fractions.
PCT_COLS = {"mean_seq_recovery", "seq_recovery", "mean_recovery", "recovery"}


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fmt(col, value):
    """Format one cell: percentages, sensible float precision, or verbatim."""
    if value is None or value == "":
        return "—"
    try:
        x = float(value)
    except ValueError:
        return str(value)
    if col in PCT_COLS:
        return f"{x * 100:.1f}%"
    if x == int(x) and abs(x) < 1e6:
        return str(int(x))
    # Trim trailing zeros so index-like columns read "0.1", not "0.100".
    return f"{x:.3f}".rstrip("0").rstrip(".")


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


def pivot_by_model(rows, index_col, value_cols):
    """Long -> wide: one row per index value, one column group per model.

    The paper's comparison tables are all "X vs. backbone-noise / temperature /
    mask fraction, full-backbone beside CA-only", which is this shape.
    """
    models = list(OrderedDict.fromkeys(r["model"] for r in rows))
    index_values = list(OrderedDict.fromkeys(r[index_col] for r in rows))
    lookup = {(r["model"], r[index_col]): r for r in rows}

    headers = [index_col]
    for m in models:
        for c in value_cols:
            label = c.replace("mean_", "").replace("_", " ")
            headers.append(f"{MODEL_LABEL.get(m, m)} {label}")

    body = []
    for iv in index_values:
        row = [str(iv)]
        for m in models:
            src = lookup.get((m, iv))
            for c in value_cols:
                row.append(fmt(c, src.get(c) if src else None))
        body.append(row)
    return md_table(headers, body)


def flat(rows, cols=None):
    """Render rows as-is, optionally restricted/ordered by `cols`."""
    cols = cols or list(rows[0].keys())
    body = [[fmt(c, r.get(c)) for c in cols] for r in rows]
    return md_table([c.replace("_", " ") for c in cols], body)


# (section title, csv path, renderer) — order matches the paper's results flow.
TABLES = [
    ("Table 1 — Temperature sweep: recovery, NLL, perplexity",
     "reproduction/phase2_reproduction/temperature_sweep_summary.csv",
     lambda rows: pivot_by_model(
         rows, "temperature", ["mean_seq_recovery", "mean_score_nll", "mean_perplexity"])),

    ("Table 2 — Input completeness: clean vs. stress-test subset",
     "reproduction/phase2_reproduction/completeness_comparison.csv",
     lambda rows: flat(rows, ["model", "temperature", "subset", "n_proteins",
                              "mean_seq_recovery", "mean_score_nll", "mean_perplexity"])),

    ("Table 3 — Robustness to Gaussian backbone noise",
     "reproduction/phase3_extensions/noise_sweep/summary.csv",
     lambda rows: pivot_by_model(
         rows, "noise", ["mean_seq_recovery", "mean_perplexity"])),

    ("Table 4 — Robustness to synthetic residue masking",
     "reproduction/phase3_extensions/masking/summary.csv",
     lambda rows: pivot_by_model(
         rows, "mask_frac_pct", ["mean_seq_recovery", "mean_perplexity"])),

    ("Table 5 — Recovery by protein length",
     "reproduction/phase3_extensions/generalization/length_bucket_summary.csv",
     lambda rows: pivot_by_model(rows, "length_bucket", ["mean_seq_recovery"])),

    ("Table 6 — Recovery by oligomeric state",
     "reproduction/phase3_extensions/generalization/monomer_multimer_summary.csv",
     lambda rows: pivot_by_model(rows, "monomer_or_multimer", ["mean_seq_recovery"])),

    ("Table 7 — Recovery by residue burial (CA-CA contact-number proxy)",
     "reproduction/phase3_extensions/generalization/buried_exposed_summary.csv",
     lambda rows: pivot_by_model(rows, "burial_class", ["mean_recovery"])),

    ("Table 8 — Low-resource inference: runtime levers",
     "reproduction/phase3_extensions/low_resource/runtime_summary.csv",
     lambda rows: flat(rows, ["lever", "num_seq_per_target", "batch_size",
                              "total_inference_s", "mean_s_per_protein", "seqs_per_min"])),

    ("Table 9 — Peak memory",
     "reproduction/phase3_extensions/low_resource/memory_summary.csv",
     lambda rows: flat(rows, ["model", "structure", "batch_size",
                              "max_rss_mb", "peak_footprint_mb"])),

    ("Table 10 — Deliberate failure cases",
     "reproduction/phase3_extensions/failure_cases/summary.csv",
     lambda rows: flat(rows, ["test", "stage", "exit_code", "mode",
                              "user_visible_warning", "detail"])),

    ("Table 11 — Per-protein runtime (Phase 2, both models, 24 seqs each)",
     "reproduction/phase2_reproduction/runtime_metrics.csv",
     lambda rows: pivot_by_model(rows, "pdb_id", ["runtime_sec", "seqs_per_sec"])),

    ("Table 12 — Dataset summary",
     "reproduction/phase1_data_prep/dataset_summary.csv",
     lambda rows: flat(rows)),
]

parts = [
    "# Paper tables",
    "",
    "Auto-generated by `reproduction/make_paper_tables.py` from the result CSVs.",
    "Do not edit by hand — re-run the script instead. Each table names the CSV",
    "it came from so any number can be traced back to the run that produced it.",
    "",
]

missing = []
for title, path, render in TABLES:
    rows = read_csv(path)
    if not rows:
        missing.append(path)
        parts += [f"## {title}", "", f"_Missing input: `{path}` — not generated._", ""]
        continue
    parts += [f"## {title}", "", f"Source: `{path}`", "", render(rows), ""]

with open(OUT_PATH, "w") as f:
    f.write("\n".join(parts) + "\n")

print(f"Wrote {len(TABLES) - len(missing)}/{len(TABLES)} tables -> {OUT_PATH}")
if missing:
    print("Missing inputs (run the corresponding phase first):")
    for p in missing:
        print("   ", p)

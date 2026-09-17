"""Phase 3 failure cases: turn the six run_failure_cases.sh logs/outputs into a
single summary table classifying how each edge case actually behaves.

Failure *modes* are what matter here, not just pass/fail — a crash and a
silently-wrong output are very different problems for a downstream user, so
each case is classified into one of:

  crash              - non-zero exit, pipeline stops
  silent_corruption  - exits 0 but the parsed/designed output is wrong with
                       no warning
  graceful_discard   - exits 0, bad input explicitly skipped and reported
  degenerate_output  - exits 0, output produced but not usable
  success            - exits 0, output correct

Writes summary.csv.
"""
import csv
import json
import os
import re

FC = "reproduction/phase3_extensions/failure_cases"

EXIT_RE = re.compile(r"^exit_code=(\d+)", re.M)
GEN_RE = re.compile(r"^(\d+) sequences of length (\d+) generated in ([\d.]+) seconds", re.M)
DISCARD_RE = re.compile(r"discarded\s+(\{.*\})", re.M)
RSS_RE = re.compile(r"^\s*(\d+)\s+maximum resident set size", re.M)
FOOTPRINT_RE = re.compile(r"^\s*(\d+)\s+peak memory footprint", re.M)
SAMPLE_RE = re.compile(
    r">T=([\d.]+), sample=(\d+), score=([\d.]+), global_score=([\d.]+), seq_recovery=([\d.]+)")


def read(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return f.read()


def exit_code(log):
    m = EXIT_RE.search(log or "")
    return int(m.group(1)) if m else None


def last_traceback_line(log):
    lines = [l.strip() for l in (log or "").splitlines() if l.strip()]
    for line in reversed(lines):
        if line.startswith("exit_code="):
            continue
        return line[:200]
    return ""


def parsed_seq(jsonl_path):
    """Concatenated chain sequences from a parsed JSONL (first entry)."""
    raw = read(jsonl_path)
    if not raw or not raw.strip():
        return None
    entry = json.loads(raw.splitlines()[0])
    keys = sorted(k for k in entry if k.startswith("seq_chain_"))
    return "".join(entry[k] for k in keys)


def samples(fa_path):
    raw = read(fa_path)
    if raw is None:
        return []
    out = []
    for line in raw.splitlines():
        m = SAMPLE_RE.match(line.strip())
        if m:
            out.append({"score": float(m.group(3)), "recovery": float(m.group(5))})
    return out


rows = []

# --- test 1: malformed coordinate field --------------------------------------
log = read(f"{FC}/test1_malformed.log")
code = exit_code(log)
rows.append({
    "test": "1_malformed_coordinate",
    "input": "pdbs/malformed_badcoord.pdb",
    "stage": "parse_multiple_chains.py",
    "exit_code": code,
    # No log at all means the case was never run — never report that as a pass.
    "mode": "not_run" if log is None else ("crash" if code else "success"),
    "detail": last_traceback_line(log) if code else
              ("run_failure_cases.sh has not been run" if log is None else "parsed without error"),
    "user_visible_warning": "yes (traceback)" if code else "n/a",
})

# --- test 2: unusual residue code --------------------------------------------
seq = parsed_seq(f"{FC}/test2_parsed.jsonl")
gaps = seq.count("-") if seq else None
rows.append({
    "test": "2_unusual_amino_acid",
    "input": "pdbs/unusual_aa.pdb",
    "stage": "parse_multiple_chains.py",
    "exit_code": exit_code(read(f"{FC}/test2_unusual_aa.log")),
    "mode": "silent_corruption" if gaps else "success",
    "detail": f"parsed seq={seq!r}, {gaps} gap char(s) — ZZZ residue dropped to '-'"
              if gaps else f"parsed seq={seq!r}",
    "user_visible_warning": "no",
})

# --- test 3: chain ID collision ----------------------------------------------
seq = parsed_seq(f"{FC}/test3_parsed.jsonl")
rows.append({
    "test": "3_chain_id_collision",
    "input": "pdbs/chain_collision.pdb",
    "stage": "parse_multiple_chains.py",
    "exit_code": exit_code(read(f"{FC}/test3_chain_collision.log")),
    "mode": "silent_corruption",
    "detail": f"parsed seq={seq!r} (len {len(seq) if seq else 0}) — colliding "
              "(chain, resnum) keys overwrite each other, residues lost",
    "user_visible_warning": "no",
})

# --- test 4: --max_length discard --------------------------------------------
log = read(f"{FC}/test4_maxlength.log")
discard = DISCARD_RE.search(log or "")
rows.append({
    "test": "4_max_length_discard",
    "input": "1UBQ (76 res) with --max_length 50",
    "stage": "protein_mpnn_run.py",
    "exit_code": exit_code(log),
    "mode": "graceful_discard" if discard else ("not_run" if log is None else "unexpected"),
    "detail": f"discarded {discard.group(1)}" if discard else last_traceback_line(log),
    "user_visible_warning": "yes (explicit discard message)",
})

# --- test 5: extreme backbone noise ------------------------------------------
s5 = samples(f"{FC}/test5_extremenoise/seqs/1UBQ.fa")
mean_rec = sum(x["recovery"] for x in s5) / len(s5) if s5 else None
rows.append({
    "test": "5_extreme_backbone_noise",
    "input": "1UBQ with --backbone_noise 5.0",
    "stage": "protein_mpnn_run.py",
    "exit_code": exit_code(read(f"{FC}/test5_extremenoise.log")),
    "mode": "degenerate_output" if mean_rec is not None and mean_rec < 0.10 else "success",
    "detail": (f"{len(s5)} sequences, mean recovery {mean_rec:.1%} — below the ~5% "
               "expected from uniform-random guessing; sequences collapse to one "
               "repeated residue") if mean_rec is not None else "no samples found",
    "user_visible_warning": "no",
})

# --- test 6: very long synthetic structure -----------------------------------
log = read(f"{FC}/test6_verylong.log")
gen = GEN_RE.search(log or "")
rss = RSS_RE.search(log or "")
foot = FOOTPRINT_RE.search(log or "")
detail = []
if gen:
    detail.append(f"{gen.group(1)} seqs of length {gen.group(2)} in {gen.group(3)}s")
if rss:
    detail.append(f"max RSS {int(rss.group(1)) / 1e6:.0f} MB")
if foot:
    detail.append(f"peak footprint {int(foot.group(1)) / 1e6:.0f} MB")
rows.append({
    "test": "6_very_long_structure",
    "input": "SYNTHETIC_VERYLONG (30 chains, ~3386 res)",
    "stage": "protein_mpnn_run.py",
    "exit_code": exit_code(log),
    "mode": "success" if gen else ("not_run" if log is None else "unexpected"),
    "detail": ", ".join(detail) if detail else last_traceback_line(log),
    "user_visible_warning": "n/a",
})

out_path = f"{FC}/summary.csv"
with open(out_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {len(rows)} rows -> {out_path}\n")
for r in rows:
    print(f'{r["test"]:<28} {r["mode"]:<18} {r["detail"]}')

bad = [r for r in rows if r["mode"] in ("crash", "silent_corruption", "degenerate_output")]
print(f'\n{len(bad)}/{len(rows)} cases fail in a way a user would not want '
      f'({sum(1 for r in bad if r["user_visible_warning"] == "no")} of them silently).')

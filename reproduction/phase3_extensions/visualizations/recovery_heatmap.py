"""Phase 3 visualization: per-residue native-vs-designed identity heatmap
for a handful of representative proteins (full-backbone model, T=0.1,
sample 1), plus a genome-browser-style match/mismatch track."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

EXAMPLE_PROTEINS = ["1UBQ", "1VII", "2GB1", "1CRN", "1ENH"]

native_seqs = {}
with open("reproduction/phase1_data_prep/native_sequences.fasta") as f:
    lines = f.read().strip().split("\n")
for i in range(0, len(lines), 2):
    pdb_id = lines[i][1:].split()[0]
    native_seqs[pdb_id] = lines[i + 1]


def load_designed(pdb_id, sample_idx=0):
    path = f"reproduction/phase2_reproduction/fullbackbone/seqs/{pdb_id}.fa"
    lines = open(path).read().strip().split("\n")
    return lines[2 * (sample_idx + 1) + 1]


fig, axes = plt.subplots(len(EXAMPLE_PROTEINS), 1, figsize=(10, 1.1 * len(EXAMPLE_PROTEINS)))
MATCH_COLOR = "#2a78d6"
MISMATCH_COLOR = "#e3e2dd"
cmap = mcolors.ListedColormap([MISMATCH_COLOR, MATCH_COLOR])

for ax, pdb_id in zip(axes, EXAMPLE_PROTEINS):
    native = native_seqs[pdb_id]
    designed = load_designed(pdb_id)
    match = np.array([[1 if n == d else 0 for n, d in zip(native, designed)]])
    ax.imshow(match, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    recovery = match.mean()
    ax.set_yticks([])
    ax.set_xlim(-0.5, len(native) - 0.5)
    ax.set_ylabel(f"{pdb_id}\n{recovery:.0%}", rotation=0, ha="right", va="center", fontsize=10)
    if pdb_id != EXAMPLE_PROTEINS[-1]:
        ax.set_xticks([])
    else:
        ax.set_xlabel("Residue position")
    for spine in ax.spines.values():
        spine.set_visible(False)

fig.suptitle("Native vs. designed sequence identity (full-backbone, T=0.1, sample 1)", y=1.02)
handles = [plt.Rectangle((0, 0), 1, 1, color=MATCH_COLOR), plt.Rectangle((0, 0), 1, 1, color=MISMATCH_COLOR)]
fig.legend(handles, ["Match", "Mismatch"], loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=2, frameon=False)
fig.tight_layout()
fig.savefig("reproduction/phase3_extensions/visualizations/recovery_heatmap.png", dpi=200, bbox_inches="tight")
print("Wrote recovery_heatmap.png")

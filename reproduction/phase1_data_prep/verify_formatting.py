"""Phase 1 verification (re-check): residue numbering continuity, insertion
codes, altloc records, and whether 'missing residues' are internal chain
breaks vs. terminal truncation, for all 30 curated PDBs."""
import glob
import os
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

warnings.simplefilter("ignore", PDBConstructionWarning)

parser = PDBParser(QUIET=True)

for path in sorted(glob.glob("inputs/phase1_dataset/pdbs/*.pdb")):
    pdb_id = os.path.basename(path).replace(".pdb", "")
    structure = parser.get_structure(pdb_id, path)
    model = next(iter(structure))

    has_insertion_code = False
    has_altloc = False
    internal_gaps = []

    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[16] != " ":
                has_altloc = True
            if line.startswith("ATOM") and line[26] != " ":
                has_insertion_code = True

    for chain in model:
        resnums = [res.id[1] for res in chain if res.id[0] == " "]
        if len(resnums) < 2:
            continue
        for i in range(1, len(resnums)):
            gap = resnums[i] - resnums[i - 1]
            if gap > 1:
                internal_gaps.append((chain.id, resnums[i - 1], resnums[i], gap - 1))

    flags = []
    if has_altloc:
        flags.append("ALTLOC")
    if has_insertion_code:
        flags.append("INSCODE")
    if internal_gaps:
        flags.append(f"INTERNAL_GAPS={internal_gaps}")

    print(f"{pdb_id}: {'OK (no formatting issues)' if not flags else ' | '.join(flags)}")

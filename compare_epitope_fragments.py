import pandas as pd
import subprocess
import re
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select

# ============================================================
# PATHS
# ============================================================

MSA_FILE = Path("aligned_cps/aligned_gene_wzg.fasta")

FRAGMENT_METADATA = Path(
    "wzg_epitope_conservation/"
    "epitope_fragments/"
    "wzg_epitope_fragment_metadata.csv"
)

UNIQUE_FRAGMENTS = Path(
    "wzg_epitope_conservation/"
    "epitope_fragments/"
    "wzg_unique_epitope_fragments.csv"
)

# Your 42 newly downloaded structures
FRAGMENT_STRUCTURE_DIR = Path(
    "wzg_epitope_conservation/"
    "epitope_fragments/"
    "predicted_structures"
)

# Change this if your original 64 Wzg structures live elsewhere
FULL_STRUCTURE_DIR = Path("pdb_outputs/pdb_outputs")

OUTPUT_DIR = Path(
    "wzg_epitope_conservation/"
    "fragment_structure_comparison"
)

EXTRACTED_DIR = OUTPUT_DIR / "native_fragments"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

TMALIGN = "TMalign"


# ============================================================
# LOAD DATA
# ============================================================

metadata = pd.read_csv(FRAGMENT_METADATA)
unique = pd.read_csv(UNIQUE_FRAGMENTS)

msa_records = list(SeqIO.parse(MSA_FILE, "fasta"))

print(f"MSA sequences: {len(msa_records)}")
print(f"Metadata rows: {len(metadata)}")
print(f"Unique fragments: {len(unique)}")


# ============================================================
# BUILD SPC -> wzg_seqN MAPPING
#
# IMPORTANT:
# The order here must correspond to the order originally used
# to generate wzg_seq0 ... wzg_seq63.
# ============================================================

spc_to_seq = {}

for i, record in enumerate(msa_records):

    spc_to_seq[record.id] = f"wzg_seq{i}"

print(f"Sequence mappings: {len(spc_to_seq)}")


# ============================================================
# RESIDUE EXTRACTION
# ============================================================

class ResidueRange(Select):

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def accept_residue(self, residue):

        residue_number = residue.id[1]

        return self.start <= residue_number <= self.end


def extract_native_fragment(
    input_cif,
    output_pdb,
    start,
    end
):
    """
    Extract a residue range from the original full-length
    Wzg structure.
    """

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "wzg",
        str(input_cif)
    )

    io = PDBIO()
    io.set_structure(structure)

    io.save(
        str(output_pdb),
        ResidueRange(start, end)
    )


# ============================================================
# TM-ALIGN
# ============================================================

def run_tmalign(native_structure, predicted_fragment):

    result = subprocess.run(
        [
            TMALIGN,
            str(native_structure),
            str(predicted_fragment)
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:

        print(result.stderr)

        raise RuntimeError(
            f"TM-align failed for {predicted_fragment}"
        )

    output = result.stdout

    # Example:
    #
    # Aligned length= 49, RMSD= 1.23,
    # Seq_ID=n_identical/n_aligned= 0.980

    match = re.search(
        r"Aligned length=\s*(\d+),"
        r"\s*RMSD=\s*([\d.]+),"
        r"\s*Seq_ID=n_identical/n_aligned=\s*([\d.]+)",
        output
    )

    if not match:
        raise ValueError(
            "Could not find alignment statistics in TM-align output."
        )

    aligned_length = int(match.group(1))
    rmsd = float(match.group(2))
    sequence_identity = float(match.group(3))

    # TM-align normally reports two scores because it normalises
    # against each structure independently.

    tm_scores = re.findall(
        r"TM-score=\s*([\d.]+)",
        output
    )

    if len(tm_scores) < 2:
        raise ValueError(
            "Could not find both TM-scores."
        )

    tm_score_native = float(tm_scores[0])
    tm_score_fragment = float(tm_scores[1])

    return {
        "aligned_length": aligned_length,
        "rmsd": rmsd,
        "sequence_identity": sequence_identity,
        "tm_score_native_normalised": tm_score_native,
        "tm_score_fragment_normalised": tm_score_fragment
    }


# ============================================================
# FIND FULL-LENGTH STRUCTURE
# ============================================================

def get_full_structure(seq_name):
    expected = FULL_STRUCTURE_DIR / f"{seq_name}.pdb"

    if expected.exists():
        return expected

    return None

# ============================================================
# RUN ALL 42 COMPARISONS
# ============================================================

results = []

for _, row in unique.iterrows():

    msa_start = int(row["msa_start"])
    msa_end = int(row["msa_end"])
    flank = int(row["flank"])
    variant = int(row["variant"])

    fragment_sequence = row["fragment_sequence"]
    representative = row["example_protein"]

    print("\n" + "=" * 60)

    print(
        f"Region {msa_start}-{msa_end} | "
        f"flank {flank} | "
        f"variant {variant}"
    )

    print(f"Representative: {representative}")

    # --------------------------------------------------------
    # Find corresponding metadata row
    # --------------------------------------------------------

    match = metadata[
        (metadata["protein"] == representative)
        & (metadata["msa_start"] == msa_start)
        & (metadata["msa_end"] == msa_end)
        & (metadata["flank"] == flank)
        & (metadata["fragment_sequence"] == fragment_sequence)
    ]

    if len(match) == 0:

        print("WARNING: No matching metadata row.")
        continue

    match = match.iloc[0]

    fragment_start = int(match["fragment_start"])
    fragment_end = int(match["fragment_end"])

    # --------------------------------------------------------
    # Find wzg_seqN
    # --------------------------------------------------------

    if representative not in spc_to_seq:

        print(
            f"WARNING: {representative} not found in MSA mapping."
        )

        continue

    seq_name = spc_to_seq[representative]

    print(f"Mapped to: {seq_name}")

    # --------------------------------------------------------
    # Full-length structure
    # --------------------------------------------------------

    full_structure = get_full_structure(seq_name)

    if full_structure is None:

        print(
            f"WARNING: Full structure missing for {seq_name}"
        )

        continue

    # --------------------------------------------------------
    # Fragment prediction
    # --------------------------------------------------------

    fragment_name = (
        f"wzg_{msa_start}_{msa_end}"
        f"_flank{flank}"
        f"_variant{variant}"
    )

    predicted_fragment = (
        FRAGMENT_STRUCTURE_DIR
        / f"{fragment_name}_model_0.cif"
    )

    if not predicted_fragment.exists():

        print(
            f"WARNING: Fragment structure missing:\n"
            f"{predicted_fragment}"
        )

        continue

    # --------------------------------------------------------
    # Extract equivalent native fragment
    # --------------------------------------------------------

    native_fragment = (
        EXTRACTED_DIR
        / f"{fragment_name}_{seq_name}_native.pdb"
    )

    print(
        f"Native residues: "
        f"{fragment_start}-{fragment_end}"
    )

    extract_native_fragment(
        full_structure,
        native_fragment,
        fragment_start,
        fragment_end
    )

    # --------------------------------------------------------
    # TM-align
    # --------------------------------------------------------

    tm = run_tmalign(
        native_fragment,
        predicted_fragment
    )

    print(
        f"TM-score: "
        f"{tm['tm_score_fragment_normalised']:.3f}"
    )

    print(
        f"RMSD: "
        f"{tm['rmsd']:.3f} Å"
    )

    results.append({

        "region": f"{msa_start}-{msa_end}",

        "msa_start": msa_start,
        "msa_end": msa_end,

        "flank": flank,
        "variant": variant,

        "representative_protein": representative,
        "wzg_sequence": seq_name,

        "fragment_start": fragment_start,
        "fragment_end": fragment_end,

        "fragment_length": len(fragment_sequence),

        **tm
    })


# ============================================================
# SAVE
# ============================================================

results_df = pd.DataFrame(results)

OUTPUT_FILE = (
    OUTPUT_DIR /
    "wzg_fragment_tm_align_results.csv"
)

results_df.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\n" + "=" * 60)

print(
    f"Completed {len(results_df)} / "
    f"{len(unique)} comparisons."
)

print(f"\nSaved to:\n{OUTPUT_FILE}")


# ============================================================
# SUMMARY
# ============================================================

if not results_df.empty:

    print("\nSummary by candidate region/flank:\n")

    summary = (
        results_df
        .groupby(["region", "flank"])
        .agg(
            n=("variant", "count"),

            mean_tm=(
                "tm_score_fragment_normalised",
                "mean"
            ),

            median_tm=(
                "tm_score_fragment_normalised",
                "median"
            ),

            mean_rmsd=(
                "rmsd",
                "mean"
            ),

            median_rmsd=(
                "rmsd",
                "median"
            )
        )
        .round(3)
    )

    print(summary)
from pathlib import Path
from Bio import SeqIO
from itertools import combinations
from config import get_gene

import pandas as pd

GENE = get_gene()

# ============================================================
# PATHS
# ============================================================

ALIGNED_DIR = Path("aligned_cps")

ALL_PROTEINS_FILE = (
    ALIGNED_DIR / "aligned_cps_proteins.fasta"
)

GENE_FILE = (
    ALIGNED_DIR / f"aligned_gene_{GENE}.fasta"
)

OUTPUT_FILE = Path(
    "complementary_gene_ranking.csv"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_serotype(record_id):
    """
    Extract the serotype from an ID such as:

        SPC06A_0004_wzg

    -> SPC06A
    """

    return record_id.split("_")[0]


def get_gene_from_filename(path):
    """
    aligned_gene_wzh.fasta -> wzh
    aligned_gene_rmlB.fasta -> rmlB
    """

    prefix = "aligned_gene_"

    return path.stem[len(prefix):]


def pairwise_identity(seq1, seq2):
    """
    Calculate pairwise amino-acid identity from two aligned sequences.

    Alignment columns where either sequence contains a gap are ignored.

    Returns a value between 0 and 1.
    """

    if len(seq1) != len(seq2):
        raise ValueError(
            "Sequences must have equal length because "
            "they should already be aligned."
        )

    matches = 0
    compared = 0

    for aa1, aa2 in zip(seq1, seq2):

        # Ignore alignment positions containing a gap
        if aa1 == "-" or aa2 == "-":
            continue

        # Ignore unknown residues
        if aa1 == "X" or aa2 == "X":
            continue

        compared += 1

        if aa1 == aa2:
            matches += 1

    if compared == 0:
        return None

    return matches / compared


# ============================================================
# FIND ALL SEROTYPES
# ============================================================

all_records = list(
    SeqIO.parse(ALL_PROTEINS_FILE, "fasta")
)

all_serotypes = {
    get_serotype(record.id)
    for record in all_records
}

print(
    f"Total serotypes in dataset: "
    f"{len(all_serotypes)}"
)


# ============================================================
# FIND SEROTYPES WITH GENE
# ============================================================

gene_records = list(
    SeqIO.parse(GENE_FILE, "fasta")
)

gene_serotypes = {
    get_serotype(record.id)
    for record in gene_records
}


# ============================================================
# FIND SEROTYPES WITHOUT GENE
# ============================================================

missing_gene_serotypes = (
    all_serotypes - gene_serotypes
)

print(
    f"Serotypes with {GENE}: "
    f"{len(gene_serotypes)}"
)

print(
    f"Serotypes without {GENE}: "
    f"{len(missing_gene_serotypes)}"
)

print(f"\n{GENE}-missing serotypes:")

print(
    ", ".join(
        sorted(missing_gene_serotypes)
    )
)


# ============================================================
# FIND ALL GENE ALIGNMENTS
# ============================================================

gene_files = sorted(
    ALIGNED_DIR.glob("aligned_gene_*.fasta")
)

print(
    f"\nGene alignments found: "
    f"{len(gene_files)}"
)


# ============================================================
# ANALYSE EACH GENE
# ============================================================

results = []


for gene_file in gene_files:

    gene = get_gene_from_filename(
        gene_file
    )

    # No point ranking Gene against itself
    if gene.lower() == GENE.lower():
        continue

    records = list(
        SeqIO.parse(gene_file, "fasta")
    )

    if not records:
        continue


    # --------------------------------------------------------
    # SEROTYPES CONTAINING THIS GENE
    # --------------------------------------------------------

    gene_serotypes = {
        get_serotype(record.id)
        for record in records
    }


    # --------------------------------------------------------
    # COVERAGE OF GENE-MISSING SEROTYPES
    # --------------------------------------------------------

    covered_missing = (
        gene_serotypes
        & missing_gene_serotypes
    )

    n_missing_covered = len(
        covered_missing
    )

    missing_coverage_percent = (
        n_missing_covered
        / len(missing_gene_serotypes)
        * 100
    )


    # --------------------------------------------------------
    # TOTAL DATASET COVERAGE
    # --------------------------------------------------------

    n_total_covered = len(
        gene_serotypes
    )

    total_coverage_percent = (
        n_total_covered
        / len(all_serotypes)
        * 100
    )


    # --------------------------------------------------------
    # GET SEQUENCES FROM ONLY THE GENE-MISSING SEROTYPES
    # --------------------------------------------------------

    missing_records = [
        record
        for record in records
        if get_serotype(record.id)
        in missing_gene_serotypes
    ]


    # --------------------------------------------------------
    # PAIRWISE SEQUENCE IDENTITY
    # --------------------------------------------------------

    identities = []

    for record1, record2 in combinations(
        missing_records,
        2
    ):

        identity = pairwise_identity(
            str(record1.seq),
            str(record2.seq)
        )

        if identity is not None:
            identities.append(identity)


    if identities:

        mean_identity = (
            sum(identities)
            / len(identities)
        )

        min_identity = min(
            identities
        )

        max_identity = max(
            identities
        )

    else:

        mean_identity = None
        min_identity = None
        max_identity = None


    # --------------------------------------------------------
    # SAVE RESULT
    # --------------------------------------------------------

    results.append({

        "gene":
            gene,

        "missing_serotypes_covered":
            n_missing_covered,

        "missing_serotypes_total":
            len(missing_gene_serotypes),

        "missing_coverage_percent":
            missing_coverage_percent,

        "total_serotypes_covered":
            n_total_covered,

        "total_serotypes":
            len(all_serotypes),

        "total_coverage_percent":
            total_coverage_percent,

        "mean_pairwise_identity":
            mean_identity,

        "min_pairwise_identity":
            min_identity,

        "max_pairwise_identity":
            max_identity,

        "n_pairwise_comparisons":
            len(identities),

        "covered_missing_serotypes":
            ",".join(
                sorted(covered_missing)
            )
    })


# ============================================================
# CREATE DATAFRAME
# ============================================================

results_df = pd.DataFrame(
    results
)


# ============================================================
# CONVERT IDENTITIES TO PERCENTAGES
# ============================================================

for column in [
    "mean_pairwise_identity",
    "min_pairwise_identity",
    "max_pairwise_identity"
]:

    results_df[column] = (
        results_df[column] * 100
    )


# ============================================================
# SORT
# ============================================================

# Primary criterion:
#   cover as many gene-missing serotypes as possible
#
# Secondary criterion:
#   sequences should be as similar as possible

results_df = results_df.sort_values(
    by=[
        "missing_serotypes_covered",
        "mean_pairwise_identity"
    ],
    ascending=[
        False,
        False
    ]
)


# ============================================================
# SAVE
# ============================================================

results_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# DISPLAY RESULTS
# ============================================================

print("\n" + "=" * 80)

print(
    "TOP COMPLEMENTARY GENES"
)

print("=" * 80)


display_columns = [

    "gene",

    "missing_serotypes_covered",

    "missing_coverage_percent",

    "total_serotypes_covered",

    "mean_pairwise_identity",

    "min_pairwise_identity"
]


print(
    results_df[
        display_columns
    ]
    .head(20)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.1f}"
    )
)


print(
    f"\nFull ranking saved to:\n"
    f"{OUTPUT_FILE}"
)

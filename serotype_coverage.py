from pathlib import Path

import pandas as pd
from Bio import SeqIO

from config import get_gene


# ============================================================
# SETTINGS
# ============================================================

GENE = get_gene()


# ============================================================
# PATHS
# ============================================================

# Contains proteins from the complete serotype dataset.
# We use this to determine the complete set of serotypes,
# including serotypes where the current gene is absent.
ALL_PROTEINS_FILE = Path(
    "aligned_cps/aligned_cps_proteins.fasta"
)

# MSA containing only proteins belonging to the current gene.
GENE_MSA_FILE = Path(
    f"aligned_cps/aligned_gene_{GENE}.fasta"
)

# Candidate epitope regions discovered by epitope_analysis.py.
CANDIDATE_FILE = Path(
    f"{GENE}_epitope_conservation/"
    f"{GENE}_candidate_regions.csv"
)

# Aggregate motif frequencies produced by motif_coverage.py.
MOTIF_FILE = Path(
    f"{GENE}_epitope_conservation/"
    f"{GENE}_motif_coverage.csv"
)

# Output from this script.
OUTPUT_FILE = Path(
    f"{GENE}_epitope_conservation/"
    f"{GENE}_serotype_coverage.csv"
)


# ============================================================
# LOAD COMPLETE SEROTYPE DATASET
# ============================================================

all_records = list(
    SeqIO.parse(
        ALL_PROTEINS_FILE,
        "fasta"
    )
)

if not all_records:
    raise ValueError(
        f"No proteins found in {ALL_PROTEINS_FILE}"
    )


# Protein names are things like:
#
# SPC01_0005_wzg
# SPC02_0005_wzh
#
# Therefore the part before the first "_" identifies
# the serotype.

all_serotypes = sorted({
    record.id.split("_")[0]
    for record in all_records
})


print("\n" + "=" * 60)
print("COMPLETE SEROTYPE DATASET")
print("=" * 60)

print(
    f"Total serotypes: "
    f"{len(all_serotypes)}"
)


# ============================================================
# LOAD GENE-SPECIFIC DATA
# ============================================================

gene_records = list(
    SeqIO.parse(
        GENE_MSA_FILE,
        "fasta"
    )
)

if not gene_records:
    raise ValueError(
        f"No {GENE} sequences found in "
        f"{GENE_MSA_FILE}"
    )


candidates = pd.read_csv(
    CANDIDATE_FILE
)

if candidates.empty:
    raise ValueError(
        f"No candidate regions found in "
        f"{CANDIDATE_FILE}"
    )


motif_summary = pd.read_csv(
    MOTIF_FILE
)

if motif_summary.empty:
    raise ValueError(
        f"No motif results found in "
        f"{MOTIF_FILE}"
    )


print(
    f"{GENE} proteins in gene MSA: "
    f"{len(gene_records)}"
)


# ============================================================
# BUILD SEROTYPE COVERAGE TABLE
# ============================================================

all_results = []


for _, candidate in candidates.iterrows():

    start = int(
        candidate["start"]
    )

    end = int(
        candidate["end"]
    )


    # --------------------------------------------------------
    # FIND DOMINANT MOTIF FOR THIS CANDIDATE
    # --------------------------------------------------------

    region_summary = motif_summary[
        (motif_summary["msa_start"] == start)
        &
        (motif_summary["msa_end"] == end)
    ]


    if region_summary.empty:

        raise ValueError(
            f"No motif summary found for "
            f"{GENE} candidate "
            f"{start}-{end}"
        )


    dominant_motif = (
        region_summary[
            "dominant_motif"
        ]
        .iloc[0]
    )


    # --------------------------------------------------------
    # MAP SEROTYPE -> GENE PROTEIN + MOTIF
    # --------------------------------------------------------

    gene_by_serotype = {}


    for record in gene_records:

        serotype = (
            record.id
            .split("_")[0]
        )

        aligned_sequence = str(
            record.seq
        )


        # Candidate coordinates are 1-indexed
        # inclusive MSA coordinates.
        #
        # For example:
        #
        # 192-200
        #
        # corresponds to:
        #
        # aligned_sequence[191:200]
        #
        # which extracts all 9 MSA positions.

        motif = aligned_sequence[
            start - 1:end
        ]


        gene_by_serotype[
            serotype
        ] = {
            "protein": record.id,
            "motif": motif
        }


    # --------------------------------------------------------
    # EXAMINE EVERY SEROTYPE
    # --------------------------------------------------------

    for serotype in all_serotypes:

        gene_data = (
            gene_by_serotype.get(
                serotype
            )
        )


        # ====================================================
        # GENE ABSENT
        # ====================================================

        if gene_data is None:

            all_results.append({
                "serotype": serotype,
                "gene": GENE,
                "gene_present": False,
                "protein": None,

                "msa_start": start,
                "msa_end": end,

                "motif": None,
                "dominant_motif":
                    dominant_motif,

                "exact_dominant_motif":
                    False,

                "motif_distance": None,

                "coverage_status":
                    "Gene absent"
            })

            continue


        # ====================================================
        # GENE PRESENT
        # ====================================================

        motif = gene_data[
            "motif"
        ]


        # ----------------------------------------------------
        # EXACT DOMINANT MOTIF?
        # ----------------------------------------------------

        exact_match = (
            motif
            == dominant_motif
        )


        # ----------------------------------------------------
        # HAMMING DISTANCE
        # ----------------------------------------------------

        # Hamming distance only makes sense here if
        # the two aligned motifs have equal lengths.

        if (
            len(motif)
            == len(dominant_motif)
        ):

            motif_distance = sum(
                aa1 != aa2
                for aa1, aa2
                in zip(
                    motif,
                    dominant_motif
                )
            )

        else:

            motif_distance = None


        # ----------------------------------------------------
        # CLASSIFY MOTIF
        # ----------------------------------------------------

        if exact_match:

            status = (
                "Exact dominant motif"
            )

        elif motif_distance == 1:

            status = (
                "Single-residue variant"
            )

        else:

            status = (
                "Other motif variant"
            )


        # ----------------------------------------------------
        # SAVE RESULT
        # ----------------------------------------------------

        all_results.append({
            "serotype": serotype,
            "gene": GENE,
            "gene_present": True,

            "protein":
                gene_data["protein"],

            "msa_start": start,
            "msa_end": end,

            "motif": motif,

            "dominant_motif":
                dominant_motif,

            "exact_dominant_motif":
                exact_match,

            "motif_distance":
                motif_distance,

            "coverage_status":
                status
        })


# ============================================================
# CREATE DATAFRAME
# ============================================================

coverage = pd.DataFrame(
    all_results
)


if coverage.empty:

    raise ValueError(
        f"No serotype coverage results "
        f"generated for {GENE}"
    )


# ============================================================
# SAVE RESULTS
# ============================================================

coverage.to_csv(
    OUTPUT_FILE,
    index=False
)


print(
    f"\nSaved serotype coverage to:\n"
    f"{OUTPUT_FILE}"
)


# ============================================================
# PRINT SUMMARY FOR EACH CANDIDATE REGION
# ============================================================

for (
    start,
    end
), region_df in coverage.groupby(
    [
        "msa_start",
        "msa_end"
    ]
):


    dominant_motif = (
        region_df[
            "dominant_motif"
        ]
        .iloc[0]
    )


    # --------------------------------------------------------
    # COUNTS
    # --------------------------------------------------------

    total_serotypes = (
        region_df[
            "serotype"
        ]
        .nunique()
    )


    gene_present = (
        region_df.loc[
            region_df[
                "gene_present"
            ],
            "serotype"
        ]
        .nunique()
    )


    gene_absent = (
        region_df.loc[
            ~region_df[
                "gene_present"
            ],
            "serotype"
        ]
        .nunique()
    )


    exact = (
        region_df.loc[
            region_df[
                "exact_dominant_motif"
            ],
            "serotype"
        ]
        .nunique()
    )


    single_variant = (
        region_df.loc[
            region_df[
                "coverage_status"
            ]
            ==
            "Single-residue variant",
            "serotype"
        ]
        .nunique()
    )


    other_variant = (
        region_df.loc[
            region_df[
                "coverage_status"
            ]
            ==
            "Other motif variant",
            "serotype"
        ]
        .nunique()
    )


    # ========================================================
    # SUMMARY
    # ========================================================

    print(
        "\n"
        + "=" * 60
    )

    print(
        f"{GENE} SEROTYPE COVERAGE: "
        f"{start}-{end}"
    )

    print(
        "=" * 60
    )


    print(
        f"\nDominant motif: "
        f"{dominant_motif}"
    )


    print(
        f"\nTotal serotypes: "
        f"{total_serotypes}"
    )


    print(
        f"{GENE} present: "
        f"{gene_present}/"
        f"{total_serotypes}"
    )


    print(
        f"{GENE} absent: "
        f"{gene_absent}/"
        f"{total_serotypes}"
    )


    print(
        f"\nExact dominant motif: "
        f"{exact}/"
        f"{total_serotypes}"
    )


    print(
        f"Single-residue variants: "
        f"{single_variant}/"
        f"{total_serotypes}"
    )


    print(
        f"Other motif variants: "
        f"{other_variant}/"
        f"{total_serotypes}"
    )


    # ========================================================
    # SEROTYPES WITHOUT GENE
    # ========================================================

    print(
        "\n"
        + "=" * 60
    )

    print(
        f"SEROTYPES WITHOUT "
        f"{GENE}"
    )

    print(
        "=" * 60
    )


    absent = region_df[
        ~region_df[
            "gene_present"
        ]
    ]


    if absent.empty:

        print("None")

    else:

        absent_serotypes = (
            absent[
                "serotype"
            ]
            .drop_duplicates()
            .tolist()
        )

        print(
            ", ".join(
                absent_serotypes
            )
        )


    # ========================================================
    # MOTIF VARIANTS
    # ========================================================

    print(
        "\n"
        + "=" * 60
    )

    print(
        f"{GENE} MOTIF VARIANTS"
    )

    print(
        "=" * 60
    )


    variants = region_df[
        region_df[
            "gene_present"
        ]
        &
        ~region_df[
            "exact_dominant_motif"
        ]
    ]


    if variants.empty:

        print("None")

    else:

        print(
            variants[
                [
                    "serotype",
                    "protein",
                    "motif",
                    "motif_distance",
                    "coverage_status"
                ]
            ]
            .to_string(
                index=False
            )
        )


    # ========================================================
    # EXACT MOTIF SEROTYPES
    # ========================================================

    print(
        "\n"
        + "=" * 60
    )

    print(
        f"SEROTYPES WITH EXACT "
        f"{dominant_motif}"
    )

    print(
        "=" * 60
    )


    exact_serotypes = (
        region_df.loc[
            region_df[
                "exact_dominant_motif"
            ],
            "serotype"
        ]
        .drop_duplicates()
        .tolist()
    )


    if exact_serotypes:

        print(
            ", ".join(
                exact_serotypes
            )
        )

    else:

        print("None")


# ============================================================
# FINISHED
# ============================================================

print(
    "\n"
    + "=" * 60
)

print(
    f"{GENE} SEROTYPE COVERAGE "
    f"ANALYSIS COMPLETE"
)

print(
    "=" * 60
)

print(
    f"\nSaved results to:\n"
    f"{OUTPUT_FILE}"
)
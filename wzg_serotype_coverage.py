import pandas as pd

from Bio import SeqIO
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

ALL_PROTEINS_FILE = Path(
    "aligned_cps/aligned_cps_proteins.fasta"
)

WZG_METADATA_FILE = Path(
    "wzg_epitope_conservation/"
    "epitope_fragments/"
    "wzg_epitope_fragment_metadata.csv"
)

OUTPUT_FILE = Path(
    "wzg_epitope_conservation/"
    "wzg_serotype_coverage.csv"
)


# ============================================================
# SETTINGS
# ============================================================

MSA_START = 370
MSA_END = 376

# Dominant motif found in our previous analysis
DOMINANT_MOTIF = "LADGDR"


# ============================================================
# GET ALL SEROTYPES IN ORIGINAL DATASET
# ============================================================

records = list(
    SeqIO.parse(ALL_PROTEINS_FILE, "fasta")
)

# Protein IDs look like:
#
# SPC01_0005_wzh
# SPC02_0006_wzd
# SPC06A_0004_wzg
#
# Therefore the first part gives us the serotype.

all_serotypes = sorted({
    record.id.split("_")[0]
    for record in records
})

print(
    f"Total serotypes in original dataset: "
    f"{len(all_serotypes)}"
)


# ============================================================
# LOAD WZG EPITOPE DATA
# ============================================================

metadata = pd.read_csv(WZG_METADATA_FILE)

region = metadata[
    (metadata["msa_start"] == MSA_START)
    & (metadata["msa_end"] == MSA_END)
][
    [
        "protein",
        "epitope_sequence"
    ]
].drop_duplicates()


print(
    f"Wzg proteins containing candidate region: "
    f"{len(region)}"
)


# ============================================================
# EXTRACT SEROTYPE FROM EACH WZG PROTEIN
# ============================================================

region["serotype"] = (
    region["protein"]
    .str.split("_")
    .str[0]
)


# ============================================================
# CHECK FOR DUPLICATE WZG PROTEINS PER SEROTYPE
# ============================================================

# Usually we expect one Wzg protein per serotype.
# This check makes sure we don't silently lose information.

duplicates = (
    region
    .groupby("serotype")
    .size()
)

duplicates = duplicates[
    duplicates > 1
]

if not duplicates.empty:

    print(
        "\nWARNING: Some serotypes have multiple "
        "Wzg proteins:"
    )

    print(duplicates)


# ============================================================
# BUILD SEROTYPE COVERAGE TABLE
# ============================================================

results = []


for serotype in all_serotypes:

    serotype_wzg = region[
        region["serotype"] == serotype
    ]

    # --------------------------------------------------------
    # NO WZG
    # --------------------------------------------------------

    if serotype_wzg.empty:

        results.append({
            "serotype": serotype,
            "wzg_present": False,
            "wzg_protein": None,
            "motif": None,
            "exact_dominant_motif": False,
            "motif_distance": None,
            "coverage_status": "No Wzg"
        })

        continue


    # --------------------------------------------------------
    # WZG PRESENT
    # --------------------------------------------------------

    for _, row in serotype_wzg.iterrows():

        motif = row["epitope_sequence"]

        exact_match = (
            motif == DOMINANT_MOTIF
        )


        # ----------------------------------------------------
        # HAMMING DISTANCE FROM DOMINANT MOTIF
        # ----------------------------------------------------

        if len(motif) == len(DOMINANT_MOTIF):

            motif_distance = sum(
                a != b
                for a, b in zip(
                    motif,
                    DOMINANT_MOTIF
                )
            )

        else:

            # Different length means ordinary Hamming
            # distance isn't defined.
            motif_distance = None


        # ----------------------------------------------------
        # STATUS
        # ----------------------------------------------------

        if exact_match:

            status = "Exact dominant motif"

        elif motif_distance == 1:

            status = "Single-residue variant"

        else:

            status = "Other motif variant"


        results.append({
            "serotype": serotype,
            "wzg_present": True,
            "wzg_protein": row["protein"],
            "motif": motif,
            "exact_dominant_motif": exact_match,
            "motif_distance": motif_distance,
            "coverage_status": status
        })


# ============================================================
# DATAFRAME
# ============================================================

coverage = pd.DataFrame(results)


# ============================================================
# SAVE
# ============================================================

coverage.to_csv(
    OUTPUT_FILE,
    index=False
)

print(
    f"\nSaved coverage table to:\n"
    f"{OUTPUT_FILE}"
)


# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("WZG SEROTYPE COVERAGE SUMMARY")
print("=" * 60)

n_serotypes = coverage["serotype"].nunique()

n_wzg = coverage.loc[
    coverage["wzg_present"],
    "serotype"
].nunique()

n_no_wzg = coverage.loc[
    ~coverage["wzg_present"],
    "serotype"
].nunique()

n_exact = coverage.loc[
    coverage["exact_dominant_motif"],
    "serotype"
].nunique()

n_single_variant = coverage.loc[
    coverage["coverage_status"]
    == "Single-residue variant",
    "serotype"
].nunique()


print(f"\nTotal serotypes: {n_serotypes}")

print(
    f"Wzg present: "
    f"{n_wzg}/{n_serotypes}"
)

print(
    f"Wzg absent: "
    f"{n_no_wzg}/{n_serotypes}"
)

print(
    f"Exact {DOMINANT_MOTIF} motif: "
    f"{n_exact}/{n_serotypes} serotypes"
)

print(
    f"Single-residue motif variant: "
    f"{n_single_variant}/{n_serotypes} serotypes"
)


# ============================================================
# SHOW SEROTYPES WITHOUT WZG
# ============================================================

print("\n" + "=" * 60)
print("SEROTYPES WITHOUT WZG")
print("=" * 60)

no_wzg = coverage[
    ~coverage["wzg_present"]
]

if no_wzg.empty:

    print("None")

else:

    for serotype in no_wzg["serotype"].unique():
        print(serotype)


# ============================================================
# SHOW NON-EXACT WZG MOTIFS
# ============================================================

print("\n" + "=" * 60)
print("WZG MOTIF VARIANTS")
print("=" * 60)

variants = coverage[
    coverage["wzg_present"]
    & ~coverage["exact_dominant_motif"]
]

if variants.empty:

    print("None")

else:

    print(
        variants[
            [
                "serotype",
                "wzg_protein",
                "motif",
                "motif_distance",
                "coverage_status"
            ]
        ].to_string(index=False)
    )


# ============================================================
# EXACT MOTIF SEROTYPES
# ============================================================

print("\n" + "=" * 60)
print(f"SEROTYPES WITH EXACT {DOMINANT_MOTIF}")
print("=" * 60)

exact = coverage[
    coverage["exact_dominant_motif"]
]

print(
    ", ".join(
        exact["serotype"].tolist()
    )
)
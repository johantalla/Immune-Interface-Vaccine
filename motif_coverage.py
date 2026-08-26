import pandas as pd
from pathlib import Path
from Bio import SeqIO

from config import get_gene


GENE = get_gene()


# ============================================================
# PATHS
# ============================================================

MSA_FILE = Path(
    f"aligned_cps/aligned_gene_{GENE}.fasta"
)

CANDIDATE_FILE = Path(
    f"{GENE}_epitope_conservation/"
    f"{GENE}_candidate_regions.csv"
)

OUTPUT_FILE = Path(
    f"{GENE}_epitope_conservation/"
    f"{GENE}_motif_coverage.csv"
)


# ============================================================
# LOAD DATA
# ============================================================

candidates = pd.read_csv(CANDIDATE_FILE)

records = list(
    SeqIO.parse(MSA_FILE, "fasta")
)

if not records:
    raise ValueError(
        f"No sequences found in {MSA_FILE}"
    )


# ============================================================
# ANALYSE EACH CANDIDATE REGION
# ============================================================

all_results = []


for _, candidate in candidates.iterrows():

    start = int(candidate["start"])
    end = int(candidate["end"])

    expected_length = end - start + 1

    print(f"\n{'=' * 60}")
    print(f"{GENE} candidate motif: {start}-{end}")
    print(f"{'=' * 60}")


    # --------------------------------------------------------
    # EXTRACT MOTIF FROM EVERY ALIGNED PROTEIN
    # --------------------------------------------------------

    motifs = []

    for record in records:

        sequence = str(record.seq)

        # MSA coordinates are 1-indexed.
        # Python slicing is 0-indexed and excludes the endpoint.
        #
        # e.g. MSA 192-200:
        # sequence[191:200]
        # gives 9 alignment positions.
        motif = sequence[start - 1:end]

        motifs.append({
            "protein": record.id,
            "msa_start": start,
            "msa_end": end,
            "motif_sequence": motif
        })


    motif_df = pd.DataFrame(motifs)


    # --------------------------------------------------------
    # COUNT MOTIF VARIANTS
    # --------------------------------------------------------

    total_proteins = motif_df["protein"].nunique()

    counts = (
        motif_df
        .groupby("motif_sequence")["protein"]
        .nunique()
        .sort_values(ascending=False)
    )


    # --------------------------------------------------------
    # DOMINANT MOTIF
    # --------------------------------------------------------

    dominant_motif = counts.index[0]
    dominant_count = counts.iloc[0]

    dominant_percentage = (
        dominant_count / total_proteins
    ) * 100


    print("\nMotif variants:")

    for motif, n_proteins in counts.items():

        percentage = (
            n_proteins / total_proteins
        ) * 100

        print(
            f"{motif}: "
            f"{n_proteins}/{total_proteins} "
            f"({percentage:.1f}%)"
        )

        all_results.append({
            "msa_start": start,
            "msa_end": end,
            "msa_length": expected_length,
            "motif_sequence": motif,
            "n_proteins": n_proteins,
            "total_proteins": total_proteins,
            "percentage": percentage,
            "dominant_motif": dominant_motif,
            "exact_dominant_motif": (
                motif == dominant_motif
            )
        })


    print(
        f"\nDominant motif: "
        f"{dominant_motif}"
    )

    print(
        f"Dominant motif coverage: "
        f"{dominant_count}/{total_proteins} "
        f"({dominant_percentage:.1f}%)"
    )


    # --------------------------------------------------------
    # SANITY CHECK
    # --------------------------------------------------------

    motif_lengths = (
        motif_df["motif_sequence"]
        .str.len()
        .unique()
    )

    if not all(
        length == expected_length
        for length in motif_lengths
    ):
        print(
            "\nWARNING: Some motifs do not have "
            f"the expected MSA length of "
            f"{expected_length}."
        )

        print(
            f"Observed lengths: "
            f"{sorted(motif_lengths)}"
        )


# ============================================================
# SAVE RESULTS
# ============================================================

if not all_results:

    raise ValueError(
        f"No motif results generated for {GENE}"
    )


results = pd.DataFrame(
    all_results
)

results.to_csv(
    OUTPUT_FILE,
    index=False
)



print(
    f"\nSaved results to:\n"
    f"{OUTPUT_FILE}"
)
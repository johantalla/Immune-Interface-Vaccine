import pandas as pd
from pathlib import Path


INPUT_FILE = Path(
    "wzg_epitope_conservation/"
    "epitope_fragments/"
    "wzg_epitope_fragment_metadata.csv"
)

OUTPUT_FILE = Path(
    "wzg_epitope_conservation/"
    "wzg_motif_coverage.csv"
)


df = pd.read_csv(INPUT_FILE)

all_results = []


for (start, end), group in df.groupby(["msa_start", "msa_end"]):

    print(f"\n{'=' * 60}")
    print(f"Wzg candidate motif: {start}-{end}")
    print(f"{'=' * 60}")

    # Each protein appears twice in the metadata because we generated
    # flank=10 and flank=20 fragments. We only want to count each
    # protein once when measuring motif coverage.
    group = group[
        [
            "protein",
            "msa_start",
            "msa_end",
            "epitope_sequence"
        ]
    ].drop_duplicates()

    total_proteins = group["protein"].nunique()

    counts = (
        group
        .groupby("epitope_sequence")["protein"]
        .nunique()
        .sort_values(ascending=False)
    )

    for sequence, n_proteins in counts.items():

        percentage = (
            n_proteins / total_proteins
        ) * 100

        print(
            f"{sequence}: "
            f"{n_proteins}/{total_proteins} proteins "
            f"({percentage:.1f}%)"
        )

        all_results.append({
            "msa_start": start,
            "msa_end": end,
            "motif_sequence": sequence,
            "n_proteins": n_proteins,
            "total_proteins": total_proteins,
            "percentage": percentage
        })


results = pd.DataFrame(all_results)

results.to_csv(
    OUTPUT_FILE,
    index=False
)

print(f"\nSaved results to:\n{OUTPUT_FILE}")
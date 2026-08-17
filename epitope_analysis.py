import pandas as pd
from Bio import SeqIO
import json

discotope = pd.read_csv(
    "wzg_epitope_conservation/wzg_discotope.csv"
)

alignment_file = "wzg_epitope_conservation/aligned_gene_wzg.fasta"

# Load JSON mapping
with open("id_mapping.json") as f:
    mapping = json.load(f)

# Load MSA
aligned_sequences = {
    record.id: str(record.seq)
    for record in SeqIO.parse(alignment_file, "fasta")
}

print("Sequences in MSA:", len(aligned_sequences))
print("Sequences in DiscoTope:", discotope["protein"].nunique())

mapped_rows = []

for protein, protein_df in discotope.groupby("protein"):

    # wzg_seq0_A -> wzg_seq0
    seq_name = protein.removesuffix("_A")

    # wzg_seq0 -> wzg_seq0.fasta
    mapping_key = seq_name + ".fasta"

    if mapping_key not in mapping:
        print(f"WARNING: {mapping_key} not found in mapping")
        continue

    # wzg_seq0.fasta -> SPC01_0005_wzg
    msa_name = mapping[mapping_key]

    if msa_name not in aligned_sequences:
        print(f"WARNING: {msa_name} not found in MSA")
        continue

    aligned_seq = aligned_sequences[msa_name]

    # Map original residue positions to MSA positions
    residue_to_msa = {}

    residue_number = 0

    for msa_position, amino_acid in enumerate(
        aligned_seq,
        start=1
    ):
        if amino_acid != "-":
            residue_number += 1
            residue_to_msa[residue_number] = msa_position

    # Map DiscoTope residues
    for _, row in protein_df.iterrows():

        res_id = int(row["res_id"])

        if res_id not in residue_to_msa:
            continue

        mapped_rows.append({
            "protein": protein,
            "serotype": msa_name.split("_")[0],
            "msa_name": msa_name,
            "res_id": res_id,
            "msa_position": residue_to_msa[res_id],
            "residue": row["residue"],
            "epitope": row["epitope"],
            "calibrated_score": row["calibrated_score"],
            "rsa": row["rsa"]
        })

mapped = pd.DataFrame(mapped_rows)

mapped.to_csv(
    "wzg_epitope_conservation/wzg_discotope_msa_mapped.csv",
    index=False
)

# --------------------------------
# Calculate epitope conservation
# --------------------------------

conservation = (
    mapped
    .groupby("msa_position")
    .agg(
        n_sequences=("protein", "nunique"),
        epitope_count=("epitope", "sum"),
        mean_discotope_score=("calibrated_score", "mean"),
        mean_rsa=("rsa", "mean")
    )
    .reset_index()
)

# Fraction of sequences where this position is predicted as an epitope
conservation["epitope_conservation"] = (
    conservation["epitope_count"]
    / conservation["n_sequences"]
)

conservation["epitope_conservation_percent"] = (
    conservation["epitope_conservation"] * 100
)

conservation.to_csv(
    "wzg_epitope_conservation/wzg_epitope_conservation.csv",
    index=False
)

print("\nHighest epitope conservation:")
print(
    conservation.sort_values(
        "epitope_conservation_percent",
        ascending=False
    ).head(20)
)

import matplotlib.pyplot as plt

plt.figure(figsize=(14, 5))

plt.plot(
    conservation["msa_position"],
    conservation["epitope_conservation_percent"]
)

plt.xlabel("Wzg MSA position")
plt.ylabel("Predicted epitope conservation (%)")

plt.ylim(0, 100)

plt.title(
    "Conservation of predicted Wzg epitopes across serotypes"
)

plt.tight_layout()

plt.savefig(
    "wzg_epitope_conservation/wzg_epitope_conservation.png",
    dpi=300
)

plt.show()

threshold = 50

high = conservation[
    conservation["epitope_conservation_percent"] >= threshold
]["msa_position"].tolist()

regions = []

if high:
    start = high[0]
    previous = high[0]

    for position in high[1:]:

        if position == previous + 1:
            previous = position

        else:
            regions.append((start, previous))
            start = position
            previous = position

    regions.append((start, previous))

print("Candidate regions:")

for start, end in regions:

    length = end - start + 1

    if length >= 5:
        region_conservation = conservation[
            conservation["msa_position"].between(start, end)
        ]

        region_mapped = mapped[
                               mapped["msa_position"].between(start,end)
        ]

        print(
            f"\n{start}-{end}",
            f"length={length}",
            f"mean epitope conservation="
            f"{region_conservation['epitope_conservation_percent'].mean():.1f}%"
        )



        for position in range(start, end + 1):

            position_df = region_mapped[
                region_mapped["msa_position"] == position
            ]

            residue_counts = position_df["residue"].value_counts()

            print(f"\nPosition {position}")
            print("Residues:")
            print(residue_counts)

            print(
                f"Epitope conservation: "
                f"{position_df['epitope'].mean() * 100:.1f}%"
            )

            print(
                f"Mean DiscoTope score: "
                f"{position_df['calibrated_score'].mean():.3f}"
            )

            print(
                f"Mean RSA: "
                f"{position_df['rsa'].mean():.3f}"
            )
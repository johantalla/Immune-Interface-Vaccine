import pandas as pd
from Bio import SeqIO
import json
import numpy as np
from pathlib import Path

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

# Store regions which pass our minimum length requirement
candidate_regions = []

print("Candidate regions:")

for start, end in regions:

    length = end - start + 1

    if length >= 5:
        # Save this region as a candidate
        candidate_regions.append((start, end))

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

candidate_df = pd.DataFrame(
    candidate_regions,
    columns=["start", "end"]
)

candidate_df.to_csv(
    "wzg_epitope_conservation/wzg_candidate_regions.csv",
    index=False
)

plddt_root = Path("boltz_outputs/predictions")
plddt_results = []

'''for start,end in candidate_regions:
    for protein,protein_df in mapped.groupby("protein"):

        seq_name = protein.removesuffix("_A")

        plddt_file = (
            plddt_root
            / seq_name
            / f"plddt_{seq_name}_model_0.npz"
        )

        if not plddt_file.exists():
            print(f"Missing plDDT file: {seq_name}")
            continue

        #Load the boltz per residue plDDT
        plddt = np.load(plddt_file)["plddt"]

        #Find the actual residues corresponding to MSA region

        region_df = protein_df[
            protein_df["msa_position"].between(start,end)
        ]

        residue_ids = region_df["res_id"].astype(int).tolist()

        if len(residue_ids) == 0:
            continue

        local_plddt = np.array([
            plddt[res_id-1]
            for res_id in residue_ids
        ])

        plddt_results.append({
            "protein": seq_name,
            "region": f"{start}-{end}",
            "n_residues": len(local_plddt),
            "mean_plddt": local_plddt.mean(),
            "median_plddt": np.median(local_plddt),
            "min_plddt": local_plddt.min(),
            "max_plddt": local_plddt.max()
        })

plddt_df = pd.DataFrame(plddt_results)

plddt_df.to_csv(
    "wzg_epitope_conservation/wzg_candidate_plddt.csv",
    index=False
)

print("\nLocal pLDDT summary:")

summary = (
    plddt_df
    .groupby("region")
    .agg(
        n_variants=("protein", "nunique"),
        mean_local_plddt=("mean_plddt", "mean"),
        median_local_plddt=("median_plddt", "median"),
        worst_variant_mean=("mean_plddt", "min"),
        best_variant_mean=("mean_plddt", "max")
    )
)'''

pae_results = []

for start, end in candidate_regions:

    for protein, protein_df in mapped.groupby("protein"):

        seq_name = protein.removesuffix("_A")

        pae_file = (
            plddt_root
            / seq_name
            / f"pae_{seq_name}_model_0.npz"
        )

        if not pae_file.exists():
            print(f"Missing PAE file: {seq_name}")
            continue

        pae = np.load(pae_file)["pae"]

        # Get actual residue IDs corresponding to this MSA region
        region_df = protein_df[
            protein_df["msa_position"].between(start, end)
        ]

        residue_ids = region_df["res_id"].astype(int).tolist()

        if len(residue_ids) == 0:
            continue

        # convert residue IDs to zero-based indices
        idx = np.array(residue_ids) - 1

        # square submatrix: all residue-vs-residue PAE values
        # within this candidate region
        local_pae = pae[np.ix_(idx, idx)]

        off_diagonal = local_pae[~np.eye(
            local_pae.shape[0],
            dtype=bool
        )]

        pae_results.append({
            "protein": seq_name,
            "region": f"{start}-{end}",
            "n_residues": len(residue_ids),
            "mean_pae": off_diagonal.mean(),
            "median_pae": np.median(off_diagonal),
            "max_pae": off_diagonal.max()
        })


pae_df = pd.DataFrame(pae_results)

pae_df.to_csv(
    "wzg_epitope_conservation/wzg_candidate_pae.csv",
    index=False
)

print("\nLocal PAE summary:")

pae_summary = (
    pae_df
    .groupby("region")
    .agg(
        n_variants=("protein", "nunique"),
        mean_local_pae=("mean_pae", "mean"),
        median_local_pae=("median_pae", "median"),
        worst_variant_mean=("mean_pae", "max"),
        best_variant_mean=("mean_pae", "min")
    )
)

print(pae_summary)

# ============================================================
# Sequence conservation vs epitope conservation
# ============================================================

plot_data = []

for position in sorted(mapped["msa_position"].unique()):

    position_df = mapped[
        mapped["msa_position"] == position
    ]

    # -------------------------
    # Sequence conservation
    # -------------------------

    residues = position_df["residue"].dropna()

    # Ignore gaps if they exist in mapped
    residues = residues[residues != "-"]

    if len(residues) > 0:
        residue_counts = residues.value_counts()

        sequence_conservation = (
            residue_counts.iloc[0] / len(residues) * 100
        )
    else:
        sequence_conservation = 0

    # -------------------------
    # Epitope conservation
    # -------------------------

    epitope_conservation = (
        position_df["epitope"].mean() * 100
    )

    plot_data.append({
        "msa_position": position,
        "sequence_conservation_percent": sequence_conservation,
        "epitope_conservation_percent": epitope_conservation
    })


comparison = pd.DataFrame(plot_data)


# ============================================================
# Plot whole Wzg protein
# ============================================================

fig, ax = plt.subplots(figsize=(16, 6))

ax.plot(
    comparison["msa_position"],
    comparison["sequence_conservation_percent"],
    label="Sequence conservation"
)

ax.plot(
    comparison["msa_position"],
    comparison["epitope_conservation_percent"],
    label="Epitope conservation"
)

# Your epitope threshold
ax.axhline(
    50,
    linestyle="--",
    label="Epitope threshold (50%)"
)

# Highlight candidate regions
for start, end in regions:

    if end - start + 1 >= 5:
        ax.axvspan(
            start,
            end,
            alpha=0.2
        )

ax.set_xlabel("MSA position")
ax.set_ylabel("Conservation (%)")

ax.set_ylim(0, 105)

ax.set_title(
    "Wzg sequence conservation vs predicted epitope conservation"
)

ax.legend()

plt.tight_layout()

plt.savefig(
    "wzg_epitope_conservation/"
    "sequence_vs_epitope_conservation.png",
    dpi=300
)

plt.show()

# ============================================================
# Zoom around each candidate region
# ============================================================

FLANK = 20

for start, end in regions:

    if end - start + 1 < 5:
        continue

    zoom_start = max(
        comparison["msa_position"].min(),
        start - FLANK
    )

    zoom_end = min(
        comparison["msa_position"].max(),
        end + FLANK
    )

    zoom = comparison[
        comparison["msa_position"].between(
            zoom_start,
            zoom_end
        )
    ]

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.plot(
        zoom["msa_position"],
        zoom["sequence_conservation_percent"],
        marker="o",
        label="Sequence conservation"
    )

    ax.plot(
        zoom["msa_position"],
        zoom["epitope_conservation_percent"],
        marker="o",
        label="Epitope conservation"
    )

    ax.axhline(
        50,
        linestyle="--",
        label="Epitope threshold (50%)"
    )

    ax.axvspan(
        start,
        end,
        alpha=0.2,
        label="Candidate epitope"
    )

    ax.set_xlabel("MSA position")
    ax.set_ylabel("Conservation (%)")

    ax.set_ylim(0, 105)

    ax.set_title(
        f"Wzg candidate region {start}-{end}"
    )

    ax.legend()

    plt.tight_layout()

    plt.savefig(
        "wzg_epitope_conservation/"
        f"conservation_zoom_{start}_{end}.png",
        dpi=300
    )

    plt.show()
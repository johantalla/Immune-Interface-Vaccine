import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import Align
import json
from config import get_gene

GENE = get_gene()

# Load TM-align results
df = pd.read_csv("tm_align_results.csv")
gene_df = df[df["gene"] == GENE].copy()

# Load mapping
with open("id_mapping.json") as f:
    mapping = json.load(f)

# Load all Wzg sequences from the ONE fasta file
sequences = {
    record.id: str(record.seq)
    for record in SeqIO.parse(f"genes/gene_{GENE}.fasta" , "fasta")
}

# Set up pairwise sequence aligner
aligner = Align.PairwiseAligner()
aligner.mode = "global"


def get_sequence(pdb_name):
    fasta_name = pdb_name.replace(".pdb", ".fasta")
    fasta_id = mapping[fasta_name]

    return sequences[fasta_id]


def sequence_identity(seq1, seq2):

    alignment = aligner.align(seq1, seq2)[0]

    counts = alignment.counts()

    return counts.identities / (
        counts.identities +
        counts.mismatches
    )


# Calculate identity for every structural comparison
identities = []

for _, row in gene_df.iterrows():

    seq1 = get_sequence(row["seq1"])
    seq2 = get_sequence(row["seq2"])

    identities.append(
        sequence_identity(seq1, seq2)
    )

gene_df["sequence_identity"] = identities


# Save results
gene_df.to_csv(
    f"{GENE}_sequence_identity_vs_tm_score.csv",
    index=False
)

print("Number of TM-align pairs:", len(gene_df))

print(gene_df[
    ["seq1", "seq2", "tm_score", "sequence_identity"]
].head(20))

print("\nTM-score range:")
print(gene_df["tm_score"].min(), gene_df["tm_score"].max())

print("\nSequence identity range:")
print(
    gene_df["sequence_identity"].min(),
    gene_df["sequence_identity"].max()
)

print("\nNaNs:")
print(gene_df[["tm_score", "sequence_identity"]].isna().sum())
# Plot

gene_df["sequence_divergence"] = 1 - gene_df["sequence_identity"]
gene_df["structural_divergence"] = 1 - gene_df["tm_score"]

plt.scatter(
    gene_df["sequence_divergence"],
    gene_df["structural_divergence"],
    alpha=0.5
)

plt.xlabel("Sequence divergence (1 - sequence identity)")
plt.ylabel("Structural divergence (1 - TM-score)")

plt.xlim(0, 0.4)
plt.ylim(0, 1)

plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


plt.figure(figsize=(8, 6))

plt.scatter(
    gene_df["sequence_identity"] * 100,
    gene_df["tm_score"],
    alpha=0.5
)

plt.xlabel("Pairwise sequence identity (%)")
plt.ylabel("TM-score")

plt.xlim(60, 100)
plt.ylim(0, 1)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig(
    f"{GENE}_sequence_identity_vs_tm_score.png",
    dpi=300
)

plt.show()

suspicious = gene_df[
    (gene_df["sequence_identity"] >= 0.95) &
    (gene_df["tm_score"] <= 0.4)
]

print(suspicious[
    ["seq1", "seq2", "sequence_identity", "tm_score"]
])

suspicious.to_csv(
    f"{GENE}_high_sequence_low_TM.csv",
    index=False
)

controls = gene_df[
    (gene_df["sequence_identity"] >= 0.95) &
    (gene_df["tm_score"] >= 0.6)
]
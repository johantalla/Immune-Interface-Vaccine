import pandas as pd
import matplotlib.pyplot as plt
import logomaker

from Bio import SeqIO
from pathlib import Path


MSA_FILE = "aligned_cps/aligned_gene_wzg.fasta"

CANDIDATE_FILE = (
    "wzg_epitope_conservation/wzg_candidate_regions.csv"
)

OUTPUT_DIR = Path("wzg_epitope_conservation/sequence_logos")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True) 

candidates = pd.read_csv(CANDIDATE_FILE)

print("Candidate regions:")
print(candidates)

records = list(SeqIO.parse(MSA_FILE, "fasta"))

print(f"\nSequences logaled: {len(records)}")

sequences = [str(record.seq) for record in records]

if not sequences:
    raise ValueError("No sequences found in hte MSA")

alignment_length = len(sequences[0])

print(f"ALignment length: {alignment_length}")

if not all(len(seq) == alignment_length for seq in sequences):
    raise ValueError( " Sequences have different lengths when the input must be aligned")

for _, candidate in candidates.iterrows():

    start = int(candidate["start"])
    end = int(candidate["end"])

    region_sequences = [seq[start -1 : end] for seq in sequences] #MSA Positions are 1-indexed hence the -1

    print(f"\nRegion {start}-{end}")
    print(f"Length:{end-start+1}")

    sequence_counts = pd.Series(region_sequences).value_counts()

    print("\nPeptide variants:")
    print(sequence_counts)


    counts_matrix = logomaker.alignment_to_matrix(sequences=region_sequences,
                                                  to_type ="counts",
                                                  characters_to_ignore =".-X"
                                                  )
    probability_matrix = counts_matrix.div(counts_matrix.sum(axis=1),axis=0)


    #### PLOTTING HERE
    region_length = end - start+1

    fig, ax = plt.subplots(
        figsize=(region_length * 1.5, 5)
    )

    logo = logomaker.Logo(
        probability_matrix,
        ax=ax,
        width=0.9,
        color_scheme="chemistry"
    )

    ax.set_title(f" Wzg candidate epitope: MSA Positions {start}-{end}")

    ax.set_xlabel("MSA position")
    ax.set_ylabel("Amino acid frequency")

    ax.set_ylim(0, 1)

    ax.set_xticks(range(region_length))
    ax.set_xticklabels(range(start, end + 1))

    plt.tight_layout()

    output_file = (
        OUTPUT_DIR /
        f"wzg_{start}_{end}_sequence_logo.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {output_file}")

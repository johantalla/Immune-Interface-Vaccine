from Bio import SeqIO
from pathlib import Path
import pandas as pd
from config import get_gene

GENE = get_gene()

#all the paths nad settings


MSA_FILE = f"aligned_cps/aligned_gene_{GENE}.fasta"

CANDIDATE_FILE = (f"{GENE}_epitope_conservation/{GENE}_candidate_regions.csv")

OUTPUT_DIR = Path(f"{GENE}_epitope_conservation/epitope_fragments")

OUTPUT_DIR.mkdir(parents=True,exist_ok=True)

FLANK_SIZES = [10,20]

#loading data

records = list(SeqIO.parse(MSA_FILE,"fasta"))
candidates = pd.read_csv(CANDIDATE_FILE)

print(f"Sequences_loaded: {len(records)}")
print("\n Candidate regions:")
print(candidates)

def msa_to_sequence_position(aligned_sequence,msa_position):
    """
    Convert a 1-indexed MSA position into a 1-indexed
    position in the original ungapped protein sequence.

    Returns None if the MSA position itself is a gap.
    """

    if aligned_sequence[msa_position - 1] == "-":
        return None

    sequence_position = sum(
        residue != "-"
        for residue in aligned_sequence[:msa_position] # just sum up all the residues until you reach that positiion discounting gaps
    )

    return sequence_position

fragment_metadata = []

for _,candidate in candidates.iterrows():
    msa_start = int(candidate["start"])
    msa_end = int(candidate["end"])

    for flank in FLANK_SIZES:

        output_file = (OUTPUT_DIR / f"{GENE}_{msa_start}_{msa_end}_flank{flank}.fasta")

        with open(output_file,"w") as handle:

            for record in records:

                aligned_sequence = str(record.seq)

                raw_sequence = aligned_sequence.replace("-","")

                seq_start = msa_to_sequence_position(aligned_sequence,msa_start)

                seq_end = msa_to_sequence_position(aligned_sequence,msa_end)

                if seq_start is None or seq_end is None:
                    print(f"SKipping {record.id}: gap at the candidate boundary {msa_start}-{msa_end}")

                fragment_start = max(1,seq_start-flank)
                fragment_end = min(len(raw_sequence),seq_end + flank)

                fragment = raw_sequence[fragment_start -1:fragment_end]

                epitope = raw_sequence[seq_start:seq_end]

                # writing the fastfile

                header = (
                    f"{record.id}"
                    f"|msa={msa_start}-{msa_end}"
                    f"|seq={seq_start}-{seq_end}"
                    f"|fragment={fragment_start}-{fragment_end}"
                    f"|flank={flank}"
                )

                handle.write(f">{header}\n")
                handle.write(f"{fragment}\n")

                fragment_metadata.append({
                    "protein": record.id,
                    "msa_start": msa_start,
                    "msa_end": msa_end,
                    "sequence_start": seq_start,
                    "sequence_end": seq_end,
                    "epitope_sequence": epitope,
                    "flank": flank,
                    "fragment_start": fragment_start,
                    "fragment_end": fragment_end,
                    "fragment_length": len(fragment),
                    "fragment_sequence": fragment
                })

        print(f"Saved: {output_file}")

        # save the metadata

metadata = pd.DataFrame(fragment_metadata)

metadata_file = (
    OUTPUT_DIR / f"{GENE}_epitope_fragment_metadata.csv"
)

metadata.to_csv(
    metadata_file,
    index=False
)

print(f"\nSaved metadata: {metadata_file}")

print("\nFragments generated:")
print(
    metadata.groupby(
        ["msa_start", "msa_end", "flank"]
    ).size()
)

print("\nFinished.")



# --------------------------------------------------
# Find unique fragment sequences
# --------------------------------------------------

unique_fragments = (
    metadata
    .groupby([
        "msa_start",
        "msa_end",
        "flank",
        "fragment_sequence"
    ])
    .agg(
        n_proteins=("protein", "count"),
        example_protein=("protein", "first")
    )
    .reset_index()
)

unique_fragments["variant"] = (
    unique_fragments
    .groupby(["msa_start", "msa_end", "flank"])
    .cumcount() + 1
)

unique_file = (
    OUTPUT_DIR / f"{GENE}_nique_epitope_fragments.csv"
)

unique_fragments.to_csv(
    unique_file,
    index=False
)

print("\nUnique fragments:")

print(
    unique_fragments.groupby(
        ["msa_start", "msa_end", "flank"]
    ).size()
)

print(f"\nSaved: {unique_file}")

# --------------------------------------------------
# Generate Boltz YAML inputs for unique fragments
# --------------------------------------------------

BOLTZ_INPUT_DIR = OUTPUT_DIR / "boltz_inputs"
BOLTZ_INPUT_DIR.mkdir(parents=True, exist_ok=True)

for _, row in unique_fragments.iterrows():

    start = int(row["msa_start"])
    end = int(row["msa_end"])
    flank = int(row["flank"])
    variant = int(row["variant"])

    sequence = row["fragment_sequence"]

    name = (
        f"{GENE}_{start}_{end}"
        f"_flank{flank}"
        f"_variant{variant}"
    )

    yaml_file = BOLTZ_INPUT_DIR / f"{name}.yaml"

    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {sequence}
"""

    with open(yaml_file, "w") as f:
        f.write(yaml_content)

print(
    f"\nGenerated {len(unique_fragments)} "
    f"Boltz input files in {BOLTZ_INPUT_DIR}"
)
from Bio import SeqIO
import pandas as pd
from glob import glob
import os
import json
with open("id_mapping.json", "r") as f:
    id_mapping = json.load(f) # boltz_filename -> original sequence id (e.g. "SPC01_0005_wzg")


def get_match_key(name): # the csv file and fasta file have different naming conventions, so we need to extract the common part of the name to match them
    parts = name.split("_")
    return "_".join(parts[:2])

def get_original_match_key(name):
    original_name = id_mapping.get("_".join(name.split("_")[:2]) + '.fasta')  # Get the original name from the mapping
    if original_name:
        parts = original_name.split("_")
        return "_".join(parts[:2])
    return None

scores_df = pd.read_csv("cps_all_proteins_merged.csv") # DiscoTope scores

scores_df["match_key"] = scores_df["protein"].apply(get_original_match_key) # Add a column to the scores dataframe for matching with fasta keys


# for each file in aligned_cps folder - perform below:

aligned_files = glob("aligned_cps/*.fasta")

# Build a fast lookup: (match_key, res_id) -> calibrated_score
score_lookup = {}
for _, row in scores_df.iterrows():
    score_lookup[(row["match_key"], row["res_id"])] = row["calibrated_score"]


alignment = list(SeqIO.parse("aligned_cps/aligned_gene_rmlB.fasta", "fasta"))

# Check what keys we're generating from this gene's FASTA
sample_ids = [rec.id for rec in alignment[:5]]
sample_keys = [get_match_key(rec.id) for rec in alignment[:5]]
print("Sample FASTA ids:", sample_ids)
print("Sample FASTA match_keys:", sample_keys)

# Check what match_keys actually exist in scores_df
print("Sample scores_df match_keys:", scores_df["match_key"].unique()[:5].tolist())

# Does even ONE of this gene's keys exist anywhere in scores_df?
found_any = scores_df["match_key"].isin(sample_keys).any()
print("Any of this gene's keys found in scores_df at all?", found_any)

for file in aligned_files:
    gene_name = os.path.basename(file).replace("aligned_", "").replace(".fasta", "")
    print(f"Processing gene: {gene_name} from file: {file}")

    alignment = list(SeqIO.parse(file, "fasta"))
    num_seqs = len(alignment)
    alignment_length = len(alignment[0].seq)


    residue_counters = {rec.id: 0 for rec in alignment} # Tracking each sequences residue counter (alignment pos -> original residue num) via rebuilding sequence as we scan columns
    fasta_keys = {rec.id: get_match_key(rec.id) for rec in alignment}


    results = []
    all_conservation = []

    for col in range(alignment_length):
        present_count = sum(1 for rec in alignment if rec.seq[col] != "-")
        all_conservation.append(present_count / num_seqs)


    for col in range(alignment_length):

        present_count = 0 # The amount of amino acids present in this column (not gaps)
        residues = [] # (res_id, original_residue_num, amino_acid)

        for rec in alignment:
            char = rec.seq[col]
            if char != "-":
                present_count += 1
                residue_counters[rec.id] += 1
                residues.append((rec.id, residue_counters[rec.id], char))
        
        conservation_percent = present_count / num_seqs

        if conservation_percent >= 0.8:
            scores_here = []
            for res_id, original_residue_num, amino_acid in residues:
                score = score_lookup.get((fasta_keys[res_id], original_residue_num))
                if score is not None:
                    scores_here.append(score)

            if scores_here:
                avg_score = sum(scores_here) / len(scores_here)
                results.append({
                    "alignment_column": col,
                    "conservation_pct": conservation_percent,
                    "num_present": present_count,
                    "avg_calibrated_score": avg_score
                })

    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df = results_df.sort_values(by="avg_calibrated_score", ascending=False)
        results_df.to_csv(f"conserved_positions_{gene_name}.csv", index=False)
    print(f"  -> {len(results_df)} candidate positions saved to conserved_positions_{gene_name}.csv")

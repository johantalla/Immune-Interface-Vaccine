import subprocess
import glob
import os
import re
import pandas as pd
import json 

def get_gene_from_boltz_name(file_name):
    base = re.sub(r'_seq\d+$', '', file_name)  # Get the gene_name part of the filename
    return base

def parse_tm_align_output(output):
    # Extract the TM-score using regex
    match = re.search(r'TM-score\s*=\s*([\d.]+)', output)
    print(f'match result: {match}')
    if match:
        return float(match.group(1))
    else:
        return None
    
def run_tm_align(seq1_path, seq2_path, output_dir):
    result = subprocess.run(['TMalign', seq1_path, seq2_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return parse_tm_align_output(result.stdout.decode('utf-8'))

pdf_files = glob.glob("pdb_outputs/pdb_outputs/*.pdb")

target_gene_csv = pd.read_csv("gene_summary.csv")

gene_groups = {}

for pdb_file in pdf_files:
    file_name = os.path.basename(pdb_file).replace(".pdb", "")
    gene_name = get_gene_from_boltz_name(file_name)    
    if gene_name:
        if gene_name not in gene_groups:
            gene_groups[gene_name] = []
        gene_groups[gene_name].append(pdb_file)

results = []

TARGET_GENES = []

for gene in target_gene_csv["gene_name"]:
    short_name = gene.split("_")[-1]
    if short_name in gene_groups:
        TARGET_GENES.append(short_name)
for gene in TARGET_GENES:
    files = gene_groups.get(gene, [])
    print(f'Processing gene: {gene} with {len(files)} PDB files')

    for i in range(len(files)):
        for j in range(i + 1, len(files)):
            seq1, seq2 = files[i], files[j]
            tm_score = run_tm_align(seq1, seq2, "tm_align_outputs")
            results.append({
                "gene": gene,
                "seq1": os.path.basename(seq1),
                "seq2": os.path.basename(seq2),
                "tm_score": tm_score
            })

results_df = pd.DataFrame(results)
results_df.to_csv("tm_align_results.csv", index=False)
print("TM-align results saved to tm_align_results.csv")
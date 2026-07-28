from Bio import SeqIO
from collections import defaultdict

gene_groups = defaultdict(list) # Dictionary to hold lists of records for each gene

for rec in SeqIO.parse("cps_proteins.seq","fasta"): # Parse the input FASTA file
    gene_name = rec.id.split("_")[-1]
    gene_groups[gene_name].append(rec) # Group records by gene name

for gene_name, records in gene_groups.items(): # Write each group to a separate FASTA file
    output_file = f"gene_{gene_name}.fasta"
    SeqIO.write(records, output_file, "fasta")
    print(f"Wrote {len(records)} sequences to {output_file}")

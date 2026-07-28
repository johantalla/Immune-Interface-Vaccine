from Bio import SeqIO
import glob
import os
import json

id_mapping = {}  # boltz_filename -> original sequence id (e.g. "SPC01_0005_wzg")

gene_files = glob.glob("genes/gene_*.fasta")
out_dir = "boltz_inputs"
os.makedirs(out_dir, exist_ok=True)

count = 0
for filepath in gene_files:
    gene_name = os.path.basename(filepath).replace("gene_", "").replace(".fasta", "")
    records = list(SeqIO.parse(filepath, "fasta"))

    for i, rec in enumerate(records):
        original_id = rec.id
        seq = str(rec.seq)

        # Boltz needs ONE sequence per file, with this exact header format:
        # >CHAIN_ID|ENTITY_TYPE|MSA_PATH
        # "empty" tells it to run in single-sequence mode (no MSA)
        boltz_filename = f"{gene_name}_seq{i}.fasta"
        out_path = os.path.join(out_dir, boltz_filename)

        with open(out_path, "w") as f:
            f.write(">A|protein|empty\n")
            f.write(f"{seq}\n")

        id_mapping[boltz_filename] = original_id
        count += 1

with open("id_mapping.json", "w") as f:
    json.dump(id_mapping, f, indent=2)

print(f"Wrote {count} individual Boltz-ready FASTA files to {out_dir}/")
print(f"Mapping saved to id_mapping.json")

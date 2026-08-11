import json

with open("id_mapping.json", "r") as f:
    mapping = json.load(f)


def get_serotype(pdb_name):
    fasta_name = pdb_name.replace(".pdb", ".fasta")
    original_name = mapping[fasta_name]
    return original_name.split("_")[0]
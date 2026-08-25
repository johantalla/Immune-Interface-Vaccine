import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from structure_utils import get_serotype
from convert_tm_score_to_distance import get_linkage
from scipy.cluster.hierarchy import dendrogram
from config import GENE

df = pd.read_csv("tm_align_results.csv")
gene_df = df[df["gene"] == GENE]

def get_serogroup(serotype_label):
    """e.g. 'SPC15A' -> '15', 'SPC01' -> '01'"""
    match = re.match(r"SPC(\d+)", serotype_label)
    return match.group(1) if match else None

Z, structures = get_linkage(gene_df)
matrix = pd.DataFrame(np.nan, index=structures, columns=structures)
for s in structures:
    matrix.loc[s, s] = 1.0
for _, row in gene_df.iterrows():
    matrix.loc[row["seq1"], row["seq2"]] = row["tm_score"]
    matrix.loc[row["seq2"], row["seq1"]] = row["tm_score"]

# Map each structure to its serotype label and serogroup
structure_to_serotype = {s: get_serotype(s) for s in structures}
structure_to_serogroup = {s: get_serogroup(structure_to_serotype[s]) for s in structures}

# Group structures by serogroup
serogroups = {}
for s, sg in structure_to_serogroup.items():
    if sg:
        serogroups.setdefault(sg, []).append(s)

# Plot one heatmap per serogroup (skip serogroups with only 1 member — nothing to compare)
for sg, members in serogroups.items():
    if len(members) < 2:
        continue

    sub_df = gene_df[gene_df["seq1"].isin(members) & gene_df["seq2"].isin(members)]
    if sub_df.empty:
        continue

    sub_Z, sub_structures = get_linkage(sub_df)
    sub_dendro = dendrogram(sub_Z, labels=sub_structures, no_plot=True)
    ordered_sub = sub_dendro["ivl"]

    sub_matrix = matrix.loc[ordered_sub, ordered_sub]
    sub_labels = [structure_to_serotype[s] for s in ordered_sub]

    plt.figure(figsize=(6, 6))
    plt.imshow(sub_matrix, cmap="viridis", origin="lower", vmin=0, vmax=1)
    plt.xticks(range(len(sub_labels)), sub_labels, rotation=90, fontsize=8)
    plt.yticks(range(len(sub_labels)), sub_labels, fontsize=8)
    cbar = plt.colorbar()
    cbar.set_label("TM-score")
    plt.title(f"Serogroup {sg} — {GENE} structural alignment")
    plt.tight_layout()
    plt.savefig(f"heatmap_serogroup_{sg}_{GENE}.png")
    plt.close()

print(f"Generated heatmaps for {sum(1 for m in serogroups.values() if len(m) >= 2)} serogroups")
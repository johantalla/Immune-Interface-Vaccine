import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from structure_utils import get_serotype
from convert_tm_score_to_distance import get_linkage
from scipy.cluster.hierarchy import dendrogram
from config import get_gene

GENE = get_gene()

df = pd.read_csv("tm_align_results.csv")

gene_df = df[df["gene"] == GENE]
Z, structures = get_linkage(gene_df)

dendro = dendrogram(
    Z,
    labels=structures,
    no_plot= True
)

ordered_structures = dendro["ivl"]

matrix = pd.DataFrame(
    np.nan,
    index=ordered_structures,
    columns=ordered_structures
)


# Diagonal = self similarity = 1
for s in ordered_structures:
    matrix.loc[s, s] = 1.0

# Fill actual TM-scores
for _, row in gene_df.iterrows():

    seq1 = row["seq1"]
    seq2 = row["seq2"]
    score = row["tm_score"]

    matrix.loc[seq1, seq2] = score
    matrix.loc[seq2, seq1] = score

# Convert structure names to serotype labels
serotype_labels = [
    get_serotype(s)
    for s in ordered_structures
]

plt.figure(figsize=(15, 15))

plt.imshow(
    matrix,
    cmap="viridis",
    origin="lower",
    vmin=0,
    vmax=1
)

plt.xticks(
    range(len(serotype_labels)),
    serotype_labels,
    rotation=90,
    fontsize=8
)

plt.yticks(
    range(len(serotype_labels)),
    serotype_labels,
    fontsize=8
)

cbar = plt.colorbar()
cbar.set_label("TM-score")

plt.title("Structural Alignment Heatmap for Wzg proteins")

plt.tight_layout()
plt.show()
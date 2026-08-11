import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from structure_utils import get_serotype


df = pd.read_csv("tm_align_results.csv")

# just one gene
gene_df = df[df["gene"] == "wzg"]

matrix = gene_df.pivot(
    index="seq1",
    columns="seq2",
    values="tm_score"
)

# fill opposite half of the matrix
matrix = matrix.combine_first(matrix.T)

# diagonal should be 1
for s in set(matrix.index) | set(matrix.columns):
    if s in matrix.index and s in matrix.columns:
        matrix.loc[s, s] = 1.0

plt.figure(figsize=(15, 15))

plt.imshow(
    matrix,
    cmap="viridis",
    origin="lower",
    vmin=0,
    vmax=1
)

plt.xticks(
    range(len(matrix.columns)),
    matrix.columns,
    rotation=90,
    fontsize=8
)

plt.yticks(
    range(len(matrix.index)),
    matrix.index,
    fontsize=8
)

cbar = plt.colorbar()
cbar.set_label("TM-score")

plt.title("Structural Alignment Heatmap for Wzg proteins")

plt.tight_layout()
plt.show()
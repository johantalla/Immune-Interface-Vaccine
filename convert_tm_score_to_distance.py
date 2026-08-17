from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import pandas as pd
import numpy as np

def get_linkage(gene_df):
    
    structures = sorted(
        set(gene_df["seq1"]) |
        set(gene_df["seq2"])
    )

    # Make distance matrix
    # Start with max distance
    dist = pd.DataFrame(
    np.ones((len(structures), len(structures))),
    index=structures,
    columns=structures
    )   

    for s in structures:
        dist.loc[s, s] = 0

        # Fill with TM-score distances
    for _, row in gene_df.iterrows():

        distance = 1 - row["tm_score"]

        dist.loc[row["seq1"], row["seq2"]] = distance
        dist.loc[row["seq2"], row["seq1"]] = distance


    # Convert to scipy condensed format
    condensed = squareform(dist.values)

    # Hierarchical clustering - use complete (discus this in notes)
    Z = linkage(
        condensed,
        method="complete",
    )
    return Z,structures
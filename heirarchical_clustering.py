
import pandas as pd
from scipy.cluster.hierarchy import dendrogram,fcluster
import matplotlib.pyplot as plt
import os
import numpy as np
from structure_utils import get_serotype
from convert_tm_score_to_distance import get_linkage

# Load CSV
df = pd.read_csv("tm_align_results.csv")
os.makedirs("gene_clusters",exist_ok=True)



for gene, gene_df in df.groupby("gene"):

    print(f"Processing {gene}")

    Z,structures = get_linkage(gene_df)

    # Assign clusters
    # TM-score threshold 0.7
    clusters = fcluster(
        Z,
        t=0.3,
        criterion="distance"
    )


    results = pd.DataFrame({
        "structure": structures,
        "cluster": clusters
    })

    results.to_csv(
        f"gene_clusters/{gene}_clusters.csv",
        index=False
    )

    serotypes = [get_serotype(s) for s in structures]


    # Plot dendrogram
    plt.figure(figsize=(8,5))

    dendrogram(
        Z,
        labels=serotypes,
        leaf_rotation=0,
        orientation='right'
    )

    plt.title(
        f"{gene} structural clustering"
    )

    plt.ylabel(
        "Distance (1 - TM-score)"
    )

    plt.tight_layout()

    plt.savefig(
        f"gene_clusters/{gene}_dendrogram.png",
        dpi=300
    )

    plt.close()
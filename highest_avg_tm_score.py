import pandas as pd
from config import get_gene

GENE = get_gene()

tm_df = pd.read_csv("tm_align_results.csv")
gene_df = tm_df[tm_df["gene"] == GENE]

# average TM-score per structure, across all its pairings
avg_scores = pd.concat([
    gene_df.groupby("seq1")["tm_score"].mean(),
    gene_df.groupby("seq2")["tm_score"].mean()
]).groupby(level=0).mean()

most_representative = avg_scores.idxmax()
print(most_representative, avg_scores.max())
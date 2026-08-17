import pandas as pd


summary = pd.read_csv("gene_summary.csv")
print(summary[summary["gene_name"] == "gene_wzg"])

wzg_candidates = pd.read_csv("conserved_positions/conserved_positions_gene_wzg.csv")
wzg_candidates = wzg_candidates.sort_values("avg_calibrated_score", ascending=False)
print(wzg_candidates.head(10))
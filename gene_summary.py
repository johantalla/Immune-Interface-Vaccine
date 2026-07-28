import pandas as pd
import glob
import os

files = glob.glob("conserved_positions_*.csv")

summary_rows = []

for f in files:
    gene_name = os.path.basename(f).replace("conserved_positions_", "").replace(".csv", "")
    df = pd.read_csv(f)

    if df.empty:
        summary_rows.append({
            "gene_name": gene_name,
            "num_candidates": 0,
            "mean_calibrated_score": None,
            "max_calibrated_score": None,
            "mean_conservation_pct": None,
            "top5_mean_scores": None
        })

        continue

    top5 = df.nlargest(5, 'avg_calibrated_score').head(5)

    summary_rows.append({
        "gene_name": gene_name,
        "num_candidates": len(df),
        "mean_calibrated_score": df['avg_calibrated_score'].mean(),
        "max_calibrated_score": df['avg_calibrated_score'].max(),
        "mean_conservation_pct": df['conservation_pct'].mean(),
        "top5_mean_scores": top5['avg_calibrated_score'].mean()
    })

summary_df = pd.DataFrame(summary_rows)
summary_df = summary_df.sort_values(by="top5_mean_scores", ascending=False)
summary_df.to_csv("gene_summary.csv", index=False)

print("Summary saved to gene_summary.csv")
print(summary_df.head(10))
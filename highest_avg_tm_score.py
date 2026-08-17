import pandas as pd



tm_df = pd.read_csv("tm_align_results.csv")
wzg_df = tm_df[tm_df["gene"] == "wzg"]

# average TM-score per structure, across all its pairings
avg_scores = pd.concat([
    wzg_df.groupby("seq1")["tm_score"].mean(),
    wzg_df.groupby("seq2")["tm_score"].mean()
]).groupby(level=0).mean()

most_representative = avg_scores.idxmax()
print(most_representative, avg_scores.max())
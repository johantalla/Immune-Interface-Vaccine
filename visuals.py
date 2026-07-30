import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("gene_summary.csv").head(50)  # Load the summary CSV and take the top 10 rows for visualization
df = df.sort_values(by="top5_mean_scores", ascending=False)

plt.figure(figsize=(10,6))
plt.barh(df['gene_name'], df['top5_mean_scores'], color='skyblue')

plt.xlabel('Top 5 Mean Calibrated Scores')
plt.ylabel('Gene Name')
plt.title('Top 5 Mean Calibrated Scores by Gene')
plt.gca().invert_yaxis()  # Invert y-axis to have the highest score at the top
plt.savefig('gene_rankings.png')
plt.show()
#spacing the y axis out
plt.subplots_adjust(left=0.2, right=0.8, top=0.99, bottom=0.01)
import pandas as pd
import glob
import os

csv_files = glob.glob("output/cps_proteins/pdb_only/*.csv")

dfs = []

for f in csv_files:

    df = pd.read_csv(f)

    protein_name = os.path.basename(f).replace("_discotope3.csv","")

    df["protein"] = protein_name

    dfs.append(df)

merged = pd.concat(dfs,ignore_index=True)

merged.to_csv("cps_all_proteins_merged.csv",index = False)


print(f"Merged {len(csv_files)} files, {len(merged)} total residues")
print(merged.head())
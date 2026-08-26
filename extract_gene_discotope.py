import pandas as pd
from pathlib import Path

from config import get_gene

GENE = get_gene()

INPUT_FILE = Path(
    "cps_all_proteins_merged.csv"
)

OUTPUT_DIR = Path(
    f"{GENE}_epitope_conservation"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_FILE = (
    OUTPUT_DIR /
    f"{GENE}_discotope.csv"
)

df = pd.read_csv(INPUT_FILE)

# protein names look like:
# wzg_seq0_A
# wzh_seq0_A
# wcjB_seq2_A
#
# so select rows belonging to this gene

gene_df = df[
    df["protein"].str.startswith(
        f"{GENE}_seq",
        na=False
    )
].copy()

if gene_df.empty:
    raise ValueError(
        f"No DiscoTope rows found for gene {GENE}"
    )

gene_df.to_csv(
    OUTPUT_FILE,
    index=False
)

print(
    f"Extracted {len(gene_df)} rows "
    f"from {gene_df['protein'].nunique()} "
    f"{GENE} proteins."
)

print(
    f"Saved to: {OUTPUT_FILE}"
)
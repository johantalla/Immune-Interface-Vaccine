## How to Run

This pipeline spans two environments: a **laptop/local machine** (for lightweight scripting, alignment, and analysis) and an **HPC cluster with GPU access** (for Boltz-2 structure prediction). DiscoTope can run on either, but is much faster on the HPC given the dataset size.
I've left some of the files I generated so you can run statistical analysis on my dataset for curiosity

### 1. Environment setup

**Boltz-2 environment (HPC, GPU required):**
```bash
python3 -m venv boltz-env
source boltz-env/bin/activate
pip install boltz
```

**DiscoTope-3.0 environment (Python 3.10 or 3.11 required):**
```bash
python3 -m venv discotope-env
source discotope-env/bin/activate
git clone https://github.com/Magnushhoie/DiscoTope-3.0/
cd DiscoTope-3.0
pip install -r requirements.txt
pip install .
unzip models.zip
```

**General analysis environment (laptop or HPC):**
Requires `pandas`, `biopython`, and `mafft` (install via `conda install -c bioconda mafft` or your system's package manager).

### 2. Prepare input sequences

Starting from a single multi-gene, multi-serotype FASTA file:

```bash
python split_by_gene.py          # splits into gene_<name>.fasta, one per gene
python prepare_boltz_inputs.py   # reformats sequences into Boltz-2's required
                                  # >CHAIN_ID|ENTITY_TYPE|MSA_PATH FASTA header format,
                                  # one sequence per file, and writes id_mapping.json
```

### 3. Run structure prediction (Boltz-2, on GPU)

Submit as an HPC batch job (adjust queue name and resource requests for your cluster):

```bash
source boltz-env/bin/activate
boltz predict boltz_inputs/ --out_dir boltz_outputs/ --accelerator gpu --use_msa_server
```

Convert output structures from mmCIF to PDB (required for DiscoTope):

```bash
python cif_to_pdb.py
```

### 4. Run epitope scoring (DiscoTope-3.0)

```bash
source discotope-env/bin/activate
python discotope3/main.py --pdb_dir pdb_outputs/ --out_dir discotope_outputs/
```

### 5. Merge DiscoTope results

```bash
python merge_cvs.py   # produces all_genes_merged.csv
```

### 6. Align sequences per gene (MAFFT)

```bash
for f in genes/gene_*.fasta; do
  name=$(basename "$f" .fasta)
  mafft --auto "$f" > "aligned_cps/aligned_${name}.fasta"
done
```

### 7. Run conservation + scoring analysis

```bash
python analysis.py
```
Produces one `conserved_positions_<gene>.csv` per gene.

### 8. Summarize and rank genes

```bash
python gene_summary.py
```
Produces `gene_summary_ranked.csv`, the final ranked output.

### Notes for reruns

- If re-running Boltz-2 with the same `--out_dir`, already-completed structures are skipped automatically.
- `id_mapping.json` is required for step 7 — do not regenerate it without also regenerating `discotope_outputs/`, as filenames must stay consistent between the two.
- Large jobs (Boltz-2 especially) should be run via the HPC job scheduler (e.g. `bsub` for LSF), not directly on a login node.

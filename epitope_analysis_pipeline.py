import argparse
import subprocess
import sys


parser = argparse.ArgumentParser()

parser.add_argument(
    "--gene",
    required=True,
    help="Gene to analyse, e.g. wzg, wzh, wze"
)

args = parser.parse_args()

GENE = args.gene


# Put these in the order they need to run.
# Replace the filenames below with the actual names of your scripts.
PIPELINE = [
    "extract_gene_discotope.py",
    "epitope_analysis.py",
    "sequence_logos.py",
    "extract_epitope_fragments.py",
    "motif_coverage.py",
    "serotype_coverage.py",
]


print("=" * 60)

print(f"RUNNING EPITOPE PIPELINE FOR: {GENE}")
print("=" * 60)


for script in PIPELINE:

    print("\n" + "=" * 60)
    print(f"RUNNING: {script}")
    print("=" * 60)

    command = [
        sys.executable,
        script,
        "--gene",
        GENE
    ]

    try:

        subprocess.run(
            command,
            check=True
        )

    except subprocess.CalledProcessError:

        print(
            f"\nERROR: {script} failed."
        )

        print(
            "Pipeline stopped so later stages "
            "do not run using incomplete results."
        )

        sys.exit(1)


print("\n" + "=" * 60)
print(f"PIPELINE COMPLETE FOR: {GENE}")
print("=" * 60)
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    "--gene",
    required=True,
    help="Gene to analyse, e.g. wzg, wzh, wzd"
)

args = parser.parse_args()

GENE = args.gene
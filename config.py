import argparse

def get_gene():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gene",
        required=True,
        help="Gene to analyse, e.g. wzg, wzh, wzd"
    )

    args = parser.parse_args()

    return args.gene
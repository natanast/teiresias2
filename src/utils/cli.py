import argparse


def parse_argument() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process a file containing amino acids."
    )
    parser.add_argument(
        "-i", "--input", type=str, help="Path to the input file containing amino acids"
    )
    return parser.parse_args()

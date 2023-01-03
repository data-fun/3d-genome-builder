"""Compute chromosome length from a genome FASTA file.

This script requires:
- a FASTA file containing the genome sequence,
- an output text file to store chromosomes sizes.
"""

import argparse

from Bio import SeqIO

def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments
    """
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "--fasta",
        action="store",
        type=str,
        help="Fasta file containing the sequence of the genome.",
        required=True,
    )
    required.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output text file containing chromosome sizes.",
        required=True,
    )
    # Add help.
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    # Parse command line arguments.
    ARGS = get_cli_arguments()

    print(f"Opening {ARGS.fasta}")
    with open(ARGS.fasta, "r") as fasta_file, open(ARGS.output, "w") as size_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            print(f"{record.id}: {len(record.seq)} bases")
            size_file.write(f"{record.id}\t{len(record.seq)}\n")

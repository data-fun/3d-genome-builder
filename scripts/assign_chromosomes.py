"""Annotate a PDB file containing a 3D genome structure with chromosome number.

To ease visualization and further analysis:
1 chromosome = 1 residue = 1 chain.

It requires:
- a PDB file containing the genome structure,
- a FASTA file containing the genome sequence,
- a Hi-C resolution.
"""

import argparse
import math
from pathlib import Path

from Bio import SeqIO
from biopandas.pdb import PandasPdb


def is_file(parser, file_path):
    """Check file exists.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Command line argument parser.
    file_path : str
        File path
    Returns
    -------
    str
        File path.
    """
    if not Path(file_path).is_file():
        parser.error(f"The file {file_path} does not exist")
    else:
        return file_path


def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments.
    """
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "--pdb",
        action="store",
        type=lambda name: is_file(parser, name),
        help="PDB file containing the 3D structure of the genome.",
        required=True,
    )
    required.add_argument(
        "--fasta",
        action="store",
        type=lambda name: is_file(parser, name),
        help="FASTA file containing the sequence of the genome.",
        required=True,
    )
    required.add_argument(
        "--resolution",
        action="store",
        type=int,
        help="HiC resolution from which the structure has been generated.",
        required=True,
    )
    required.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output PDB file containing the annotated 3D structure of the genome.",
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


def extract_chromosome_length(fasta_name):
    """Extract chromosome length from a FASTA file.

    Parameters
    ----------
    fasta_name : str
        Name of Fasta file containing the sequence of the genome.

    Returns
    -------
    list
        List of chromosome lengths.
    """
    chromosome_length_lst = []
    with open(fasta_name, "r") as fasta_file:
        print(f"Reading {fasta_name}")
        for record in SeqIO.parse(fasta_file, "fasta"):
            length = len(record.seq)
            print(f"Found chromosome {record.id} with {length} bases")
            chromosome_length_lst.append(length)
    return chromosome_length_lst


def assign_chromosome_number(
    pdb_name_in, chromosome_length, HiC_resolution, pdb_name_out
):
    """Assign chromosome numbers to the whole genome 3D structure.

    Note:
    - The PDB file produced by Pastis is not readable by Biopython
    because the residue number column is missing.
    - We use instead the biopandas library: http://rasbt.github.io/biopandas/

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome.
    chromosome_length : list
        List with chromosome lengths.
    HiC_resolution : int
        HiC resolution.
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome.
    """
    coordinates = PandasPdb().read_pdb(pdb_name_in)
    coordinates_df = coordinates.df["ATOM"]
    # Verify the number of beads in the structure
    # and the number of beads deduced from the sequence and resolution
    # are the same.
    print(f"Number of beads read from structure: {coordinates_df.shape[0]}")
    beads_per_chromosome = [
        math.ceil(length / HiC_resolution) for length in chromosome_length
    ]
    print(
        f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}"
    )
    if coordinates_df.shape[0] != sum(beads_per_chromosome):
        raise ValueError(
            "Number of beads in structure and "
            "number of beads deduced from sequence and resolution "
            "are different."
        )
    # Add residue (chromosome) number.
    residue_number = []
    for index, bead_number in enumerate(beads_per_chromosome):
        residue_number = residue_number + [index + 1] * bead_number
    coordinates_df["residue_number"] = residue_number

    # Add residue name based on residue number.
    coordinates_df["residue_name"] = coordinates_df["residue_number"].apply(
        lambda x: f"C{x:02d}"
    )

    # Add chain (chromosome) name.
    # Chromosome 1 -> Chain A
    # Chromosome 2 -> Chain B
    # ...
    # Chromosome 26 -> Chain Z
    # Chromosome 27 -> Chain A
    # Chromosome 28 -> Chain B
    # ...
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    # Handle case where there are more than 26 chromosomes.
    repeat_factor = math.ceil(len(chromosome_length) / len(letters))
    residue_to_chain = dict(
        zip(coordinates_df["residue_number"].unique(), letters * repeat_factor)
    )
    # Map residue number to chain id.
    coordinates_df["chain_id"] = coordinates_df["residue_number"].map(residue_to_chain)

    # Fix atom number starting at 1 (instead of 0).
    coordinates_df["atom_number"] = coordinates_df["atom_number"] + 1

    # Write new PDB file to disk.
    coordinates.df["ATOM"] = coordinates_df
    coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read FASTA file and get chromosome length.
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number.
    assign_chromosome_number(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.resolution, ARGS.output)

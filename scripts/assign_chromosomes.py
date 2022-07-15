"""Annotate a PDB file containing a 3D genome structure with chromosome number.

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
        Command line argument parser
    file_path : str
        File path
    Returns
    -------
    str
        File path    
    """
    print(type(parser))
    if not Path(file_path).is_file():
        parser.error(f"The file {file_path} does not exist")
    else:
        return file_path


def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb",
        action="store",
        type=lambda name: is_file(parser, name),
        help="PDB file containing the 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "--fasta",
        action="store",
        type=lambda name: is_file(parser, name),
        help="Fasta file containing the sequence of the genome",
        required=True,
    )
    parser.add_argument(
        "--resolution",
        action="store",
        type=int,
        help="HiC resolution",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output PDB file containing the annotated 3D structure of the genome",
        required=True,
    )
    return parser.parse_args()


def extract_chromosome_length(fasta_name):
    """Extract chromosome length from a FASTA file.

    Parameters
    ----------
    fasta_name : str
        Name of Fasta file containing the sequence of the genome
    
    Returns
    -------
    list
        List of chromosome lengthes
    """
    chromosome_length_lst = []
    with open(fasta_name, "r") as fasta_file:
        print(f"Reading {fasta_name}")
        for record in SeqIO.parse(fasta_file, "fasta"):
            length = len(record.seq)
            print(f"Found chromosome {record.id} with {length} bases")
            chromosome_length_lst.append(length)
    return chromosome_length_lst


def assign_chromosome_number(pdb_name_in, chromosome_length, HiC_resolution, pdb_name_out):
    """Assign chromosome numbers to the whole genome 3D structure.

    Note:
    - The PDB file produced by Pastis is not readable by Biopython
    because the residue number column is missing.
    - We use instead the biopandas library. http://rasbt.github.io/biopandas/

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome
    chromosome_length : list
        List with chromosome lengths
    HiC_resolution : int
        HiC resolution
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    beads_per_chromosome = [math.ceil(length/HiC_resolution) for length in chromosome_length]
    print(f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}")
    
    # Add residue number
    residue_number = []
    for index, bead_number in enumerate(beads_per_chromosome):
        residue_number = residue_number + [index+1] * bead_number
    pdb_coordinates.df["ATOM"]["residue_number"] = residue_number

    # Add residue name
    pdb_coordinates.df["ATOM"]["residue_name"] = "CHR"

    # Fix atom number starting at 1 (instead of 0)
    pdb_coordinates.df["ATOM"]["atom_number"] = pdb_coordinates.df["ATOM"]["atom_number"] + 1 

    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    assign_chromosome_number(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.resolution, ARGS.output)
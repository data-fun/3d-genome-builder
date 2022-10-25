"""Convert a PDB file containing a 3D genome structure to a g3d file.

It requires:
- a PDB file containing the genome structure,
- a FASTA file containing the genome sequence,
- a Hi-C resolution.
"""

import argparse
import math
import pandas as pd
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
        "--resolution", action="store", type=int, help="HiC resolution", required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output g3d file containing the annotated 3D structure of the genome",
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


def convert_to_g3d(pdb_name_in, chromosome_length, HiC_resolution, g3d_name_out):
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
    g3d_name_out : str
        Output g3d file containing the annotated 3D structure of the genome
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    beads_per_chromosome = [
        math.ceil(length / HiC_resolution) for length in chromosome_length
    ]
    print(
        f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}"
    )

    # Deduce chromosomal position of each atoms
    bp_coordinates = list()
    for i in range(len(chromosome_length)):
        bp_coordinates_chrom_x = list(range(0, chromosome_length[i], HiC_resolution))
        bp_coordinates = bp_coordinates + bp_coordinates_chrom_x

    # Extract chromosomes and 3d coordinates information
    g3d_data = { 'chrom': pdb_coordinates.df["ATOM"]["residue_number"],
                'locus':bp_coordinates,
                'X3D_x': pdb_coordinates.df["ATOM"]["x_coord"],
                'X3D_y': pdb_coordinates.df["ATOM"]["y_coord"],
                'X3D_z': pdb_coordinates.df["ATOM"]["z_coord"] }
    result = pd.DataFrame(g3d_data)

    #save as g3d file
    result.to_csv(g3d_name_out, sep="\t", index=False)

if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    convert_to_g3d(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.resolution, ARGS.output)
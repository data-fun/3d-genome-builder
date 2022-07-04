"""Annotate a PDB file containing a 3D genome structure with a quantitative parameter.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a resolution.
"""

import argparse
import math
import pandas as pd

from Bio import SeqIO
from biopandas.pdb import PandasPdb



def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--pdb",
        action="store",
        type=str,
        help="PDB file containing the 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--fasta",
        action="store",
        type=str,
        help="Fasta file containing the sequence of the genome",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--annotation",
        action="store",
        type=str,
        help="csvfile containing the quantitative parameter",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--resolution",
        action="store",
        type=int,
        help="HiC resolution",
        required=True,
    )
    parser.add_argument(
        "-o",
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


def map_parameter(pdb_name_in, chromosome_length, annotation, HiC_resolution, pdb_name_out):
    """Assign a quantitative parameter to the whole genome 3D structure.

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
    annotation : str
        txt file containing the annotation of a quantitative parameter
    HiC_resolution : int
        HiC resolution
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    beads_per_chromosome = [math.ceil(length/HiC_resolution) for length in chromosome_length]
    print(f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}")

    pdb_coordinates_df = pdb_coordinates.df["ATOM"]

    chipseq = pd.read_csv(annotation, sep="\t", header=None)
    chrom_names = {"chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4, "chr5": 5, "chr6": 6, "chr7": 7}
    chipseq[0]= chipseq[0].map(chrom_names)

    pdb_coordinates_df["b_factor"]=chipseq[3]

    pdb_coordinates.df["ATOM"] = pdb_coordinates_df
    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    map_parameter(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.annotation, ARGS.resolution, ARGS.output)
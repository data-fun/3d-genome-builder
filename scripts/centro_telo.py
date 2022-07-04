"""Annotate a PDB file containing a 3D genome structure with centromeres and telomeres.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a txt file containing the centromere and telomere positions,
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
        help="txt file containing the annotation of centromeres and telomeres",
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


def del_index(index, list):
    outliers_bp_coordinates_chrom_x = []
    for ele in sorted(index, reverse = True):
        outliers_bp_coordinates_chrom_x.append(list[ele])
        del list[ele]
    return [list, outliers_bp_coordinates_chrom_x]


def assign_centromere_telomere(pdb_name_in, chromosome_length, annotation, HiC_resolution, pdb_name_out):
    """Annotate centromeres and telomeres to the whole genome 3D structure : C for centromeres, T for telomeres and S for sequence.

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
        txt file containing the annotation of centromeres and telomeres
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
    pdb_coordinates_df_output = pd.DataFrame()

    # Read annotation file
    centromeres = pd.read_csv(annotation, sep="\t", header=None)

    for i in range(len(chromosome_length)):
        chrom_x = pdb_coordinates_df[pdb_coordinates_df["residue_number"]==i+1]
        chrom_x.reset_index(drop=True, inplace=True)
        chrom_x = chrom_x.assign(atoms_bp_coor=pd.Series(list(range(0, chromosome_length[i], HiC_resolution))))
        

        # Annotate beads included in centromeres and telomeres
        chrom_x["chain_id"] = "S"
        chrom_x.loc[chrom_x["atoms_bp_coor"].between(centromeres.iloc[i, 1], centromeres.iloc[i, 2]), "chain_id"] = "C"
        print(str(sum(chrom_x["chain_id"]=="C"))+" beads associated to the centromere of chromosome "+str(i+1))
        
        pdb_coordinates_df_output = pd.concat([pdb_coordinates_df_output, chrom_x.iloc[:,:-1]], axis=0)

    
    pdb_coordinates.df["ATOM"] = pdb_coordinates_df_output
    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    assign_centromere_telomere(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.annotation, ARGS.resolution, ARGS.output)
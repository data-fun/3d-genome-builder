"""Interpolate genes from a PDB file containing a 3D genome structure.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a csv file containing the genes positions,
- a resolution.
"""

import argparse
import math
import pandas as pd

from Bio import SeqIO
from biopandas.pdb import PandasPdb
from numpy import NaN
from scipy.interpolate import PchipInterpolator


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
        help="txt file containing the annotation of genes",
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


def interpolate_genes(pdb_name_in, chromosome_length, annotation, HiC_resolution, pdb_name_out):
    """Interpolate genes according to a PDB file containing a 3D genome structure.

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
        csv file containing the annotation of genes
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
    pdb_coordinates_df_output = pdb_coordinates_df.iloc[:1,]

    loci_bp_coordinates = pd.read_csv(annotation, sep="\t", header=None)
    loci_bp_coordinates = loci_bp_coordinates[loci_bp_coordinates[2]=="CDS"]
    loci_bp_coordinates[0].replace({"Supercontig_14.1": 1,
                                        "Supercontig_14.2": 2,
                                        "Supercontig_14.3": 3,
                                        "Supercontig_14.4": 4,
                                        "Supercontig_14.5": 5,
                                        "Supercontig_14.6": 6,
                                        "Supercontig_14.7": 7}, inplace=True)
    loci_bp_coordinates["gene_name"] = loci_bp_coordinates[8].str[9:17]
    loci_bp_coordinates.reset_index(drop=True, inplace=True)
    loci_bp_coordinates.drop([1, 2, 5, 6, 7, 8], axis = 1, inplace=True)
    loci_bp_coordinates.rename(columns = {0:"chrom", 3:"start", 4:"stop"}, inplace = True)

    new_atoms = pd.DataFrame(columns = ["residue_number", "x_coord", "y_coord", "z_coord"])

    for i in range(len(chromosome_length)):
        loci_bp_coordinates_chrom_x = loci_bp_coordinates[loci_bp_coordinates["chrom"]==i+1]
        loci_bp_coordinates_chrom_x = loci_bp_coordinates_chrom_x["start"]
        chrom_x = pdb_coordinates_df[pdb_coordinates_df["residue_number"]==i+1]
        chrom_x = [list(chrom_x["x_coord"]), list(chrom_x["y_coord"]), list(chrom_x["z_coord"])]

        atoms_bp_coor = list(range(0, chromosome_length[i], HiC_resolution))
        x_3d_coor = chrom_x[0]
        loci_x_3d_coor = PchipInterpolator(atoms_bp_coor, x_3d_coor)(loci_bp_coordinates_chrom_x)
        y_3d_coor = chrom_x[1]
        loci_y_3d_coor = PchipInterpolator(atoms_bp_coor, y_3d_coor)(loci_bp_coordinates_chrom_x)
        z_3d_coor = chrom_x[2]
        loci_z_3d_coor = PchipInterpolator(atoms_bp_coor, z_3d_coor)(loci_bp_coordinates_chrom_x)

        new_atoms_chrom_x = pd.DataFrame(zip([i+1]*len(loci_x_3d_coor), loci_x_3d_coor, loci_y_3d_coor, loci_z_3d_coor), columns=["residue_number", "x_coord", "y_coord", "z_coord"])
        new_atoms = pd.concat([new_atoms, new_atoms_chrom_x], sort=False)

    new_atoms.reset_index(drop=True, inplace=True)
    pdb_coordinates_df_output = pd.concat([pdb_coordinates_df_output.iloc[:-1,], new_atoms], axis=0)
    pdb_coordinates_df_output = pdb_coordinates_df_output.assign(record_name="ATOM",
                                                                    atom_number=list(range(0,len(new_atoms))),
                                                                    blank_1="",
                                                                    atom_name="O",
                                                                    alt_loc="",
                                                                    residue_name="CHR",
                                                                    blank_2="",
                                                                    chain_id="G",
                                                                    insertion="",
                                                                    blank_3="",
                                                                    occupancy=1,
                                                                    b_factor=75,
                                                                    blank_4="",
                                                                    segment_id="",
                                                                    element_symbol="",
                                                                    charge=NaN,
                                                                    line_idx=pdb_coordinates_df_output.index)

    pdb_coordinates.df["ATOM"] = pdb_coordinates_df_output
    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    interpolate_genes(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.annotation, ARGS.resolution, ARGS.output)
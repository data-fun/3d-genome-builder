"""Interpolate nan 3D coordinates in a PDB file containing a 3D genome structure.

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


def interpolate_nan_coordinates(pdb_name_in, chromosome_length, HiC_resolution, pdb_name_out):
    """Interpolate nan 3D coordinates in a PDB file containing a 3D genome structure.

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

    # Tag beads with nan coordinates as "missing"
    pdb_coordinates_df = pdb_coordinates.df["ATOM"]
    pdb_coordinates_missing = pdb_coordinates_df.assign(missing = "point")
    pdb_coordinates_missing.loc[pdb_coordinates_missing["x_coord"].isnull(), "missing"] = "missing"
    index_whole_genome = pdb_coordinates_missing.loc[pdb_coordinates_missing["missing"]=="missing"].index.values.tolist()

    for i in range(len(chromosome_length)):

        # Select index of missing values from one chromosome
        chrom_x = pdb_coordinates_missing[pdb_coordinates_missing["residue_number"]==i+1]
        chrom_x.reset_index(drop=True, inplace=True)
        atoms_bp_coor = list(range(0, chromosome_length[i], HiC_resolution))
        index = chrom_x.loc[chrom_x["missing"]=="missing"].index.values.tolist()

        if index != []:
            print(str(len(index))+" missing beads in chromosome "+str(i+1)+" !")

            # Interpolate missing beads coordinates from the other beads coordinates
            atoms_bp_coor, missing_bp_coordinates_chrom_x = del_index(index, atoms_bp_coor)
            chrom_x = chrom_x[chrom_x["missing"]=="point"]
            chrom_x_list = [list(chrom_x["x_coord"]), list(chrom_x["y_coord"]), list(chrom_x["z_coord"])]
            x_3d_coor = chrom_x_list[0]
            missing_new_x_3d_coor = PchipInterpolator(atoms_bp_coor, x_3d_coor)(missing_bp_coordinates_chrom_x)
            y_3d_coor = chrom_x_list[1]
            missing_new_y_3d_coor = PchipInterpolator(atoms_bp_coor, y_3d_coor)(missing_bp_coordinates_chrom_x)
            z_3d_coor = chrom_x_list[2]
            missing_new_z_3d_coor = PchipInterpolator(atoms_bp_coor, z_3d_coor)(missing_bp_coordinates_chrom_x)
            zipped = list(zip(missing_new_x_3d_coor, missing_new_y_3d_coor, missing_new_z_3d_coor))
            new_atoms = pd.DataFrame(zipped, columns=["x_coord", "y_coord", "z_coord"], index=index_whole_genome[:len(index)])
            index_whole_genome = index_whole_genome[len(index):]
            
            # Save interpolated coordonates
            pdb_coordinates_df["x_coord"] = pdb_coordinates_df["x_coord"].fillna(new_atoms["x_coord"])
            pdb_coordinates_df["y_coord"] = pdb_coordinates_df["y_coord"].fillna(new_atoms["y_coord"])
            pdb_coordinates_df["z_coord"] = pdb_coordinates_df["z_coord"].fillna(new_atoms["z_coord"])

    pdb_coordinates.df["ATOM"] = pdb_coordinates_df
    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)

    # Assign chromosome number
    interpolate_nan_coordinates(ARGS.pdb, CHROMOSOME_LENGTH, ARGS.resolution, ARGS.output)

"""Interpolate nan 3D coordinates in a PDB file containing a 3D genome structure.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a resolution.
"""

import argparse
import pandas as pd

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
        "-o",
        "--output",
        action="store",
        type=str,
        help="Output PDB file containing the annotated 3D structure of the genome",
        required=True,
    )
    return parser.parse_args()


def del_index(index, list):
    outliers_bp_coordinates_chrom_x = []
    for ele in sorted(index, reverse = True):
        outliers_bp_coordinates_chrom_x.append(list[ele])
        del list[ele]
    return [list, outliers_bp_coordinates_chrom_x]


def interpolate_nan_coordinates(pdb_name_in, pdb_name_out):
    """Interpolate nan 3D coordinates in a PDB file containing a 3D genome structure.

    Note:
    - The PDB file produced by Pastis is not readable by Biopython
    because the residue number column is missing.
    - We use instead the biopandas library. http://rasbt.github.io/biopandas/

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    # Tag beads with nan coordinates as "missing"
    pdb_coordinates_df = pdb_coordinates.df["ATOM"]
    pdb_coordinates_missing = pdb_coordinates_df.assign(missing = "point")
    pdb_coordinates_missing.loc[pdb_coordinates_missing["x_coord"].isnull(), "missing"] = "missing"
    index_whole_genome = pdb_coordinates_missing.loc[pdb_coordinates_missing["missing"]=="missing"].index.values.tolist()
    chromosomes_number = max(pdb_coordinates_df["residue_number"])

    for i in range(chromosomes_number):
        # Select index of missing values from one chromosome
        chrom_x = pdb_coordinates_missing[pdb_coordinates_missing["residue_number"]==i+1]
        chrom_x.reset_index(drop=True, inplace=True)
        index_known_beads = chrom_x.loc[chrom_x["missing"]=="point"].index.values.tolist()
        index = chrom_x.loc[chrom_x["missing"]=="missing"].index.values.tolist()

        if index != []:
            print(str(len(index))+" missing beads in chromosome "+str(i+1)+" !")

            # Interpolate missing beads coordinates from the other beads coordinates
            chrom_x = chrom_x[chrom_x["missing"]=="point"]
            chrom_x_list = [list(chrom_x["x_coord"]), list(chrom_x["y_coord"]), list(chrom_x["z_coord"])]
            x_3d_coor = chrom_x_list[0]
            missing_new_x_3d_coor = PchipInterpolator(index_known_beads, x_3d_coor)(index)
            y_3d_coor = chrom_x_list[1]
            missing_new_y_3d_coor = PchipInterpolator(index_known_beads, y_3d_coor)(index)
            z_3d_coor = chrom_x_list[2]
            missing_new_z_3d_coor = PchipInterpolator(index_known_beads, z_3d_coor)(index)
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

    # Assign chromosome number
    interpolate_nan_coordinates(ARGS.pdb, ARGS.output)

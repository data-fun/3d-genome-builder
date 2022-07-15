"""Interpolate missing coordinates in a PDB file containing a 3D genome structure.

This script requires:
- a PDB file containing the genome structure.
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
        "--pdb",
        action="store",
        type=str,
        help="PDB file containing the 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output PDB file containing the fixed 3D structure of the genome",
        required=True,
    )
    return parser.parse_args()


def interpolate_missing_coordinates(pdb_name_in, pdb_name_out):
    """Interpolate missing coordinates in PDB file containing a 3D genome structure.

    3D coordinates of the missing beads are interpolated one by one.

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome
    """
    pdb = PandasPdb().read_pdb(pdb_name_in)
    atoms = pdb.df["ATOM"]
    print(f"Found {atoms.shape[0]} beads in {pdb_name_in}")
    for residue_number in atoms["residue_number"].unique():
        missing_beads = (
            atoms.query(f"residue_number == {residue_number}")[
                ["x_coord", "y_coord", "z_coord"]
            ]
            .isna()
            .sum()
            .sum()
        ) // 3
        if not missing_beads:
            continue
        print(f"Chromosome {residue_number} has {missing_beads} missing beads")
        chromosome = atoms.query(f"residue_number == {residue_number}")
        beads_index = chromosome.dropna(
            subset=["x_coord", "y_coord", "z_coord"]
        ).index.to_numpy()
        beads_x = chromosome.dropna(subset=["x_coord", "y_coord", "z_coord"])[
            "x_coord"
        ].to_numpy()
        beads_y = chromosome.dropna(subset=["x_coord", "y_coord", "z_coord"])[
            "y_coord"
        ].to_numpy()
        beads_z = chromosome.dropna(subset=["x_coord", "y_coord", "z_coord"])[
            "z_coord"
        ].to_numpy()
        missing_beads_index = chromosome.loc[
            chromosome[["x_coord", "y_coord", "z_coord"]].isna().any(1)
        ].index.to_numpy()
        # Build models for each chromosome and each coordinate
        interpolate_x = PchipInterpolator(beads_index, beads_x)
        interpolate_y = PchipInterpolator(beads_index, beads_y)
        interpolate_z = PchipInterpolator(beads_index, beads_z)
        for bead_index in missing_beads_index:
            atoms.loc[bead_index, "x_coord"] = interpolate_x(bead_index)
            atoms.loc[bead_index, "y_coord"] = interpolate_y(bead_index)
            atoms.loc[bead_index, "z_coord"] = interpolate_z(bead_index)
    # atoms["atom_name"] = "CA"
    # atoms["element_symbol"] = "C"
    # atoms["b_factor"] = 0.0
    pdb.df["ATOM"] = atoms
    pdb.to_pdb(pdb_name_out)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    interpolate_missing_coordinates(ARGS.pdb, ARGS.output)

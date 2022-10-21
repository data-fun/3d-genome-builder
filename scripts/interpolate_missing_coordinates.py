"""Interpolate missing coordinates in a PDB file containing a 3D genome structure.

This script requires:
- a PDB file containing the genome structure.
"""

import argparse
import pandas as pd

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
    interpolate_atoms = pd.DataFrame(columns=atoms.columns)
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
        chromosome = atoms.query(f"residue_number == {residue_number}")
        if not missing_beads:
            interpolate_atoms = pd.concat([interpolate_atoms, chromosome])
            continue
        print(f"Chromosome {residue_number} has {missing_beads} missing beads")
        
        # Build models for each chromosome
        #TODO faire l'interpolation en une seule ligne, sur les 3 colonnes
        interpolate_chromosome = chromosome.copy()
        interpolate_chromosome["x_coord"] = chromosome["x_coord"].interpolate(method='pchip', axis=0, limit_area='inside')
        interpolate_chromosome["y_coord"] = chromosome["y_coord"].interpolate(method='pchip', axis=0, limit_area='inside')
        interpolate_chromosome["z_coord"] = chromosome["z_coord"].interpolate(method='pchip', axis=0, limit_area='inside')
        deleted_atoms = interpolate_chromosome["x_coord"].isna().sum()
        print(f"Removed {deleted_atoms} atoms from chromosome {residue_number} extremities")
            
        interpolate_atoms = pd.concat([interpolate_atoms, interpolate_chromosome])
    # atoms["atom_name"] = "CA"
    # atoms["element_symbol"] = "C"
    # atoms["b_factor"] = 0.0
    interpolate_atoms.dropna(subset=["x_coord", "y_coord", "z_coord"], inplace=True)
    print(f"Found {interpolate_atoms.shape[0]} beads after interpolation and extremities nan beads removal.")
    
    pdb.df["ATOM"] = interpolate_atoms
    pdb.to_pdb(pdb_name_out)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    interpolate_missing_coordinates(ARGS.pdb, ARGS.output)

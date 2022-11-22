"""Add missing beads in the 3D structure of the genome.

Interpolate missing coordinates in a PDB file containing a 3D genome structure.
"""

import argparse

from biopandas.pdb import PandasPdb


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
        "--input-pdb",
        action="store",
        type=str,
        help="Input PDB file containing the 3D structure of the genome.",
        required=True,
    )
    required.add_argument(
        "--output-pdb",
        action="store",
        type=str,
        help="Output PDB file containing the completed 3D structure of the genome.",
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
    # Interpolate missing coordinates chromosome by chromosome.
    coord_columns= ["x_coord", "y_coord", "z_coord"]
    for residue_number in atoms["residue_number"].unique():
        missing_beads = (
            atoms.query(f"residue_number == {residue_number}")
            [coord_columns]
            .isna()
            .sum()   # Sum NAs over columns (coordinates).
            .sum()   # Sum NAs over lines.
        ) // 3       # 1     missing atom has 3 missing coordinates.
        chromosome = atoms.query(f"residue_number == {residue_number}")
        if not missing_beads:
            continue
        print(f"Chromosome {residue_number} has {missing_beads} missing beads")

        # Interpolate coordinates on the three coordinates (x, y, z).
        # Only interpolation between existing beads.
        # No extrapolation at chromosome ends.
        interpolate_chromosome = chromosome.copy()
        interpolate_chromosome[coord_columns] = chromosome[coord_columns].interpolate(
            method="pchip", axis=0, limit_area="inside"
        )
        # Print number of missing beads that will later remove.
        deleted_atoms_number = interpolate_chromosome["x_coord"].isna().sum()
        print(
            f"Removed {deleted_atoms_number} beads from chromosome {residue_number} extremities"
        )
        # Save interpolated chromosome into the full structure.
        atoms.loc[atoms["residue_number"] == residue_number, coord_columns] = interpolate_chromosome[coord_columns]
    # Remove in one step all beads with missing coordinates
    # located at chromosome extremities.
    atoms = atoms.dropna(subset=coord_columns)
    print(
        f"Found {atoms.shape[0]} beads after interpolation and removal of extremities missing beads."
    )

    pdb.df["ATOM"] = atoms
    pdb.to_pdb(pdb_name_out)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    interpolate_missing_coordinates(ARGS.input_pdb, ARGS.output_pdb)

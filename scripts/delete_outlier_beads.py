"""Delete outliers beads in the 3D structure of the genome.

Delete outliers coordinates in a PDB file containing a 3D genome structure.
"""

import argparse
import numpy as np
import pandas as pd

from biopandas.pdb import PandasPdb

# Allow to print wide dataframes.
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)


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
    optional.add_argument(
        "--threshold",
        action="store",
        type=float,
        help="Distance threshold to detect flipped contigs. Default: 3",
        required=False,
        default=2.0,
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


def delete_outlier_beads(pdb_name_in, pdb_name_out, threshold):
    """Delete outlier coordinates in a PDB file containing a 3D genome structure.

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

    coord_columns = ["x_coord", "y_coord", "z_coord"]
    ATOMS = pd.DataFrame()

    # delete outlier coordinates chromosome by chromosome.
    for residue_number in atoms["residue_number"].unique():
        # Select beads of one chromosome
        chromosome_df = atoms.query(f"residue_number == {residue_number}").reset_index(
            drop=True
        )

        # Compute Euclidean distances between bead n and bead n+1
        chromosome_df["distance"] = np.linalg.norm(
            chromosome_df[coord_columns].to_numpy()
            - chromosome_df[coord_columns].shift(-1).to_numpy(),
            axis=1,
        )
        # The last distance is Nan because there is no bead n+1 for the last bead.
        # We replace it by the distance between bead n and bead n-1.
        column_index = chromosome_df.columns.get_loc("distance")
        chromosome_df.iat[-1, column_index] = chromosome_df.iat[-2, column_index]

        # Compute median distance with possible Nan values
        median_distance = chromosome_df["distance"].median(skipna=True)
        print(f"Median distance between beads: {median_distance:.2f}")

        # Select outlier beads.
        # i.e. the beads with the greater ditances.
        # Output beads coordinates with distances
        outlier_beads = chromosome_df["distance"] >= threshold * median_distance

        if not outlier_beads.sum():
            print(f"Chromosome {residue_number} has no outlier beads")
            chromosome_df = chromosome_df.drop("distance", axis=1)
            # Save chromosome into the full stucture.
            ATOMS = pd.concat([ATOMS, chromosome_df])

        else:
            # Print number of outlier beads that will be later removed.
            print(
                f"Chromosome {residue_number} has {outlier_beads.sum()} outlier beads"
            )
            deleted_atoms = chromosome_df[
                chromosome_df["distance"] >= threshold * median_distance
            ]
            print(
                deleted_atoms[
                    [
                        "record_name",
                        "atom_number",
                        "residue_name",
                        "residue_number",
                        "x_coord",
                        "y_coord",
                        "z_coord",
                        "distance",
                    ]
                ]
            )
            kept_atoms = chromosome_df[
                chromosome_df["distance"] < threshold * median_distance
            ]
            kept_atoms = kept_atoms.drop("distance", axis=1)
            # Save chromosome into the full stucture.
            ATOMS = pd.concat([ATOMS, kept_atoms])
    ATOMS.reset_index(inplace=True, drop=True)
    ATOMS["line_idx"] = ATOMS.index
    print(f"Found {ATOMS.shape[0]} beads after removal of outlier beads.")
    pdb.df["ATOM"] = ATOMS
    pdb.to_pdb(pdb_name_out)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    delete_outlier_beads(ARGS.input_pdb, ARGS.output_pdb, ARGS.threshold)

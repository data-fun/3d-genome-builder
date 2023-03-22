"""Annotate a PDB file containing a 3D genome structure with a quantitative parameter.

This parameter could be issued from ChIP-Seq experiment for instance.

It requires:
- a PDB file with the genome structure,
- a Bedgraph file with a quantitative parameter in the last column,
  see https://genome.ucsc.edu/goldenPath/help/bedgraph.html.
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
        help="PDB file with the 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "--bedgraph",
        action="store",
        type=str,
        help="Bedgraph file with the quantitative parameter",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output PDB file with the annotated 3D structure of the genome",
        required=True,
    )
    return parser.parse_args()


def normalize_parameter(values, max_value=999):
    """Normalize the quantitative parameter values.

    B factors values must be between 0 and 999.99.

    We normalize the values between 0 and max_value,
    only if the max value is above max_value.
    """
    if max(values) > max_value:
        values = values / max(values) * max_value
    return values


def map_parameter(pdb_name_in, bedgraph_name, pdb_name_out):
    """Assign a quantitative parameter to a genome 3D structure.

    Parameters
    ----------
    pdb_name_in : str
        PDB file with the 3D structure of the genome
    bedgraph_name : str
        Bedgraph file with a quantitative parameter
    pdb_name_out : str
        Output PDB file with the annotated 3D structure of the genome
    """
    # Read the PDB file.
    pdb = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb.df['ATOM'].shape[0]}")
    atoms = pdb.df["ATOM"]

    # Read the bedGraph file.
    quantitative_parameter = pd.read_csv(bedgraph_name, sep="\t", header=None)
    quantitative_parameter.columns = ["chromosome", "start", "end", "value"]

    # Get the index of the beads.
    # Having deleted missing is not an issue.
    atom_number = list(atoms["atom_number"] - 1)

    # Delete the quantitative parameter values corresponding to deleted beads.
    quantitative_parameter = quantitative_parameter.iloc[atom_number]
    quantitative_parameter.reset_index(drop=True, inplace=True)

    # Normalize the quantitative parameter values,
    # only is max value is above 999.
    quantitative_parameter["value"] = normalize_parameter(
        quantitative_parameter["value"]
    )

    # Map the quantitative parameter values to the corresponding beads, in the "b_factor" pdb column.
    atoms["b_factor"] = quantitative_parameter["value"]

    # Write the PDB file.
    pdb.df["ATOM"] = atoms
    pdb.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()
    map_parameter(ARGS.pdb, ARGS.bedgraph, ARGS.output)

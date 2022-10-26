"""Annotate a PDB file containing a 3D genome structure with a quantitative parameter.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a BedGraph file containning the quantitative parameter,
- a resolution.
"""

import argparse
from symbol import atom
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
        "--pdb",
        action="store",
        type=str,
        help="PDB file containing the 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "--BedGraph",
        action="store",
        type=str,
        help="BedGraph file containing the quantitative parameter",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        type=str,
        help="Output PDB file containing the annotated 3D structure of the genome",
        required=True,
    )
    return parser.parse_args()


def map_parameter(pdb_name_in, BedGraph, pdb_name_out):
    """Assign a quantitative parameter to the whole genome 3D structure.

    Note:
    - The PDB file produced by Pastis is not readable by Biopython
    because the residue number column is missing.
    - We use instead the biopandas library. http://rasbt.github.io/biopandas/

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome
    BedGraph : str
        BedGraph file containing the annotation of a quantitative parameter
    HiC_resolution : int
        HiC resolution
    pdb_name_out : str
        Output PDB file containing the annotated 3D structure of the genome
    """

    pdb = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb.df['ATOM'].shape[0]}")
    atoms = pdb.df["ATOM"]
    quantitative_parameter = pd.read_csv(BedGraph, sep="\t", header=None)

    # Get the index of the beads without the ones already removed from the pdb file (by interpolate_missing_coordinates.py )
    atom_number = list(atoms["atom_number"]-1)

    # Delete the quantitative parameter values corresponding to removed beads
    quantitative_parameter = quantitative_parameter.iloc[atom_number]
    quantitative_parameter.reset_index(drop=True, inplace=True)

    # Map the quantitative parameter values to the corresponding beads, in the "b_factor" pdb column
    #TODO ajouter la normalisation pour les valeurs au dessus de 999.99
    atoms["b_factor"]=quantitative_parameter[3]

    pdb.df["ATOM"] = atoms
    pdb.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Assign chromosome number
    map_parameter(ARGS.pdb, ARGS.BedGraph, ARGS.output)
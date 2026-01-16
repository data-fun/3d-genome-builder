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
        help="Output PDB file with the annotated 3D structure of the genome (a corresponding mmCIF file with the same base name is also written).",
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

def write_mmcif_from_pandas(atoms, cif_name):
    """Write a minimal mmCIF file from a Pandas DataFrame.
    
    The function exports atom/bead coordinates and associated B-factor values
    to a simplified _atom_site mmCIF format, intended for visualization of
    3D genome structures in tools such as Mol*.

    Parameters
    ----------
    atoms : pandas.DataFrame
        DataFrame containing atomic or bead-level coordinates and annotations.
    cif_name : str
        Path to the output mmCIF file.
    """
    out = open(cif_name, "w")
    out.write("#\nloop_\n")
    out.write("_atom_site.group_PDB\n")
    out.write("_atom_site.label_seq_id\n")
    out.write("_atom_site.label_atom_id\n")
    out.write("_atom_site.label_asym_id\n")
    out.write("_atom_site.label_comp_id\n")
    out.write("_atom_site.Cartn_x\n")
    out.write("_atom_site.Cartn_y\n")
    out.write("_atom_site.Cartn_z\n")
    out.write("_atom_site.B_iso_or_equiv\n")

    to_print = atoms.to_string(
        header=False, index=False,
        columns=[
            "record_name",
            "atom_number",
            "atom_name",    
            "chain_id",
            "residue_name",
            "x_coord",
            "y_coord",
            "z_coord",
            "b_factor"
        ]
    )

    out.write(to_print)
    out.close()

def map_parameter(pdb_name_in, bedgraph_name, pdb_name_out):
    """Assign a quantitative parameter to a genome 3D structure.

    The function assigns values from a bedGraph file (e.g. ChIP-seq signal)
    to the B-factor field of a PDB-based 3D genome model. The annotated
    structure is written both as a PDB file and as a minimal mmCIF file,
    enabling interactive browsing of genomic regions in Mol*, which is not
    supported by the PDB format alone.

    Parameters
    ----------
    pdb_name_in : str
        PDB file with the 3D structure of the genome
    bedgraph_name : str
        Bedgraph file with a quantitative parameter
    pdb_name_out : str
        Output PDB file with the annotated 3D structure of the genome
        A corresponding mmCIF file is also written.
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
    pdb.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    pdb.df["ATOM"] = atoms
    print(f"Wrote {pdb_name_out}")

    # Write the mmCIF file.
    cif_name_out = pdb_name_out.replace(".pdb", ".cif")
    write_mmcif_from_pandas(atoms, cif_name_out)
    print(f"Wrote {cif_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()
    map_parameter(ARGS.pdb, ARGS.bedgraph, ARGS.output)

"""Detect and flip inverted contigs from a PDB file containing a 3D genome structure.

Flip inverted contigs in the 3D structure and in the sequence of the genome.

This script requires:
- a PDB file containing the 3D genome structure,
- a FASTA file containing the genome sequence,
- the HiC resolution from which the structure has been generated.
"""

import argparse
import math
from operator import invert
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd


def get_cli_arguments():
    """Command line argument parser.

    Returns
    -------
    argparse.Namespace
        Object containing arguments
    """
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "--pdb",
        action="store",
        type=str,
        help="PDB file containing the 3D structure of the genome.",
        required=True,
    )
    required.add_argument(
        "--fasta",
        action="store",
        type=str,
        help="Fasta file containing the sequence of the genome.",
        required=True,
    )
    required.add_argument(
        "--resolution",
        action="store",
        type=int,
        help="HiC resolution from which the structure has been generated.",
        required=True,
    )
    required.add_argument(
        "--output-pdb",
        action="store",
        type=str,
        help="Output PDB file containing the fixed 3D structure of the genome.",
        required=True,
    )
    required.add_argument(
        "--output-fasta",
        action="store",
        type=str,
        help="Output FASTA file containing the fixed sequence of the genome.",
        required=True,
    )
    optional.add_argument(
        "--threshold",
        action="store",
        type=float,
        help="Distance threshold to detect flipped contigs. Default: 3",
        required=False,
        default=3.0,
    )
    optional.add_argument(
        "--debug",
        action="store_true",
        help="Debug flag. Output one TSV file per chromosome with bead distances.",
        required=False,
        default=False,
    )
    # Trigger to run the script or not.
    optional.add_argument(
        "--run",
        action="store",
        choices=("True", "False"),
        help="Run the script. Default: False.",
        required=False,
        default="False",
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


def extract_chromosome_name_length(fasta_name):
    """Extract chromosome name and length from a FASTA file.

    Parameters
    ----------
    fasta_name : str
        Name of the FASTA file containing the sequence of the genome.

    Returns
    -------
    tuple
        List of chromosome names.
        List of chromosome lengthes.
    """
    chromosome_name_lst = []
    chromosome_length_lst = []
    with open(fasta_name, "r") as fasta_file:
        print(f"Reading genome sequence in {fasta_name}")
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id
            length = len(record.seq)
            print(f"Found chromosome {name} with {length} bases")
            chromosome_name_lst.append(name)
            chromosome_length_lst.append(length)
    return chromosome_name_lst, chromosome_length_lst


def find_inverted_contigs(
    pdb_name_in, chromosome_lengths, HiC_resolution, threshold
):
    """Find inverted contigs.

    Inverted contigs are detected based on the eucledian distance
    between adjacent beads in the 3D structure of the genome.

    Parameters
    ----------
    pdb_name_in : str
        PDB file containing the 3D structure of the genome.
    chromosome_lengths : list
        List with chromosome lengths.
    HiC_resolution : int
        HiC resolution.
    threshold : float
        Threshold to detect flipped contigs.

    Returns
    -------
    inverted_contigs : dict
        Dictionnary with inverted contigs.
    """
    pdb_structure = PandasPdb().read_pdb(pdb_name_in)
    structure_df = pdb_structure.df["ATOM"]
    print(f"Number of beads read from structure: {structure_df.shape[0]}")

    coord_columns = ["x_coord", "y_coord", "z_coord"]
    missing_beads_number = structure_df[coord_columns].isna().sum().sum() // 3
    print(
        f"Found {missing_beads_number} missing beads in structure"
    )

    beads_per_chromosome = [
        math.ceil(length / HiC_resolution) for length in chromosome_lengths
    ]
    print(
        f"Number of expected beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}"
    )

    if structure_df.shape[0] != sum(beads_per_chromosome):
        sys.exit(
            f"Cannot process structure {pdb_name_in} because it contains "
            f"{structure_df.shape[0]} beads "
            f"instead of {sum(beads_per_chromosome)}"
        )

    inverted_contigs = {}

    for chrom_num in structure_df["residue_number"].unique():
        print(f"\nLooking for inverted contigs into chromosome {chrom_num}")

        # Select beads of one chromosome
        chromosome_df = structure_df.query(
            f"residue_number == {chrom_num}"
        ).reset_index(drop=True)

        # Compute Euclidean distances between bead n and bead n+1
        chromosome_df["distance"] = np.linalg.norm(
            chromosome_df[coord_columns].to_numpy() - chromosome_df[coord_columns].shift(-1).to_numpy(),
            axis=1
        )
        # Compute median distance with possible Nan values
        median_distance = chromosome_df["distance"].median(skipna=True)
        print(f"Median distance between beads: {median_distance:.2f}")

        # Select extremities of inverted contigs
        # i.e. beads with distance above a given threshold.
        # Output beads coordinates with distances
        if ARGS.debug:
            filename = f"chr_{chrom_num}.tsv"
            target_columns = [
                "record_name", "atom_number", "atom_name",
                "residue_name", "chain_id", "residue_number",
                "x_coord", "y_coord", "z_coord", "line_idx",
                "distance",
            ]
            print(f"DEBUG: writing {filename} with distances.")
            chromosome_df[target_columns].to_csv(
                filename, sep="\t", index=False, na_rep="nan"
            )
        
        beads_selection = (
            chromosome_df["distance"] > threshold * median_distance
        )
        inversion_limits = chromosome_df.loc[
            beads_selection, "atom_number"
        ].values
        if len(inversion_limits) % 2 != 0:
            print("WARNING: odd number of inversion limits found")
            print(
                "WARNING: this might lead to a wrong detection of inverted contigs"
            )
            print(inversion_limits)
        if len(inversion_limits) != 0:
            for limit_1, limit_2 in zip(
                inversion_limits[0::2], inversion_limits[1::2]
            ):
                print(
                    f"Chromosome {chrom_num}: found inverted contig between bead {limit_1+1} and bead {limit_2}"
                )
                if chrom_num in inverted_contigs:
                    inverted_contigs[chrom_num].append((limit_1 + 1, limit_2))
                else:
                    inverted_contigs[chrom_num] = [(limit_1 + 1, limit_2)]
        else:
            inverted_contigs[chrom_num] = []
    return inverted_contigs


def flip_inverted_contigs_in_structure(
    inverted_contigs, pdb_name_in, pdb_name_out
):
    """Flip inverted contigs in the 3D structure of the genome.

    Parameters
    ----------
    inverted_contigs : dict
        Dictionnary with inverted contigs
    pdb_name_in : str
        PDB file containing the 3D structure of the genome
    pdb_name_out : str
        Output PDB file containing the 3D structure of the genome
    """
    pdb_structure = PandasPdb().read_pdb(pdb_name_in)
    coordinates = pdb_structure.df["ATOM"]
    if sum(map(len, inverted_contigs.values()))== 0:
        pdb_structure.to_pdb(
        path=pdb_name_out, records=None, gz=False, append_newline=True
    )
        print("\nNothing to fix. Structure is fine.")
        return
    print("\nFlipping contigs.")
    for chrom_num in inverted_contigs:
        for contig in inverted_contigs[chrom_num]:
            contig_start, contig_end = contig
            print(
                f"Structure of chromosome {chrom_num}: "
                f"flip contig between beads {contig_start} "
                f"and {contig_end}"
            )
            contig_start_index = coordinates[
                (coordinates["residue_number"] == chrom_num)
                & (coordinates["atom_number"] == contig_start)
            ].index[0]
            contig_end_index = coordinates[
                (coordinates["residue_number"] == chrom_num)
                & (coordinates["atom_number"] == contig_end)
            ].index[0]
            contig_before_df = coordinates.loc[: contig_start_index - 1, :]
            contig_df = coordinates.loc[contig_start_index:contig_end_index, :]
            contig_after_df = coordinates.loc[contig_end_index + 1 :, :]
            # Flip contig.
            contig_df = contig_df[::-1]
            # Assemble genome structure.
            coordinates = pd.concat(
                [contig_before_df, contig_df, contig_after_df]
            )

    coordinates = coordinates.reset_index(drop=True)
    # The 'line_idx' column keeps the real order of atoms in the PDB file.
    coordinates["line_idx"] = coordinates.index
    pdb_structure.df["ATOM"] = coordinates
    pdb_structure.to_pdb(
        path=pdb_name_out, records=None, gz=False, append_newline=True
    )


def flip_inverted_contigs_in_sequence(
    inverted_contigs,
    chromosome_names,
    fasta_name_in,
    HiC_resolution,
    fasta_name_out,
):
    """Flip inverted contigs in the genome 3D structure and sequence.

    Parameters
    ----------
    inverted_contigs : dict
        Dictionnary with inverted contigs.
    chromosome_names : list
        List with chromosome names.
    fasta_name_in : str
        Name of Fasta file containing the sequence of the genome.
    HiC_resolution : int
        HiC resolution.
    fasta_name_out : str
        Output FASTA file containing the fixed sequence (at the 3D structure resolution!).
    """
    # Flip inverted contigs in the genome sequence.
    genome_fasta = SeqIO.to_dict(SeqIO.parse(fasta_name_in, "fasta"))

    if bool([a for a in inverted_contigs.values() if a == []]):
        with open(fasta_name_out, "w") as fasta_file:
            SeqIO.write(genome_fasta.values(), fasta_file, "fasta")
        print("Nothing to fix. Sequence is fine.")
        return

    for chrom_num in inverted_contigs:
        chrom_name = chromosome_names[chrom_num - 1]
        chrom_sequence = str(genome_fasta[chrom_name].seq)
        for contig in inverted_contigs[chrom_num]:
            contig_start = contig[0] * HiC_resolution
            contig_end = contig[1] * HiC_resolution
            print(
                f"Sequence of chromosome {chrom_num}: "
                f"flip inverted contig between base {contig_start} "
                f"and {contig_end}"
            )
            contig_sequence = chrom_sequence[contig_start : contig_end + 1]
            # Flip contig.
            contig_sequence = contig_sequence[::-1]
            # Reassemble chromosome sequence.
            chrom_sequence = (
                chrom_sequence[:contig_start]
                + contig_sequence
                + chrom_sequence[contig_end + 1 :]
            )
        genome_fasta[chrom_name].seq = Seq(chrom_sequence)

    # Write genome sequence.
    with open(fasta_name_out, "w") as fasta_file:
        SeqIO.write(genome_fasta.values(), fasta_file, "fasta")


if __name__ == "__main__":
    # Parse command line arguments.
    ARGS = get_cli_arguments()

    if ARGS.run != "True":
        print("Do not verify inverted contigs.")
        print(f"Copying {ARGS.pdb} to {ARGS.output_pdb}.")
        sys.exit()

    # Read Fasta file and extract chromosome names and lengths.
    CHROMOSOME_NAMES, CHROMOSOME_LENGTHS = extract_chromosome_name_length(
        ARGS.fasta
    )

    # Find inverted contigs.
    INVERTED_CONTIGS = find_inverted_contigs(
        ARGS.pdb, CHROMOSOME_LENGTHS, ARGS.resolution, ARGS.threshold
    )
    # Flip inverted contigs in the genome 3D structure and sequence.
    flip_inverted_contigs_in_structure(
        INVERTED_CONTIGS, ARGS.pdb, ARGS.output_pdb
    )
    flip_inverted_contigs_in_sequence(
        INVERTED_CONTIGS,
        CHROMOSOME_NAMES,
        ARGS.fasta,
        ARGS.resolution,
        ARGS.output_fasta,
    )

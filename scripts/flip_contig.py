"""Modify a PDB file containing a 3D genome structure to flip inverted contigs. The corrected sequence (at the 3D structure resolution!) is saved in a new fasta file.

It requires:
- a PDB file containing the genome structure,
- a fasta file containing the genome sequence,
- a resolution.
"""

import argparse
import math
import numpy as np
import pandas as pd

from Bio import SeqIO
from biopandas.pdb import PandasPdb
from Bio.Seq import Seq

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
        "-o_pdb",
        "--output_pdb",
        action="store",
        type=str,
        help="Output PDB file containing the annotated 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "-o_fasta",
        "--output_fasta",
        action="store",
        type=str,
        help="Output FASTA file containing the corrected sequence of the genome",
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


def extract_chromosome_name(fasta_name):
    """Extract chromosome name from a FASTA file.

    Parameters
    ----------
    fasta_name : str
        Name of Fasta file containing the sequence of the genome
    
    Returns
    -------
    list
        List of chromosome names
    """
    chromosome_name_lst = []
    with open(fasta_name, "r") as fasta_file:
        print(f"Reading {fasta_name}")
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id
            print(f"Found chromosome {name}")
            chromosome_name_lst.append(name)
    return chromosome_name_lst


def flip_inverted_contig(pdb_name_in, chromosome_length, chromosome_name, fasta_name, HiC_resolution, pdb_name_out, fasta_name_out):
    """Flip contigs inverted in the genome assembly according to the Euclidian distances between beads in the 3D genome structure.

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
    chromosome_name : list
        List with chromosome names
    fasta_name : str
        Name of Fasta file containing the sequence of the genome
    HiC_resolution : int
        HiC resolution
    pdb_name_out : str
        Output PDB file containing the 3D structure of the genome
    fasta_name_out : str
        Output FASTA file containing the corrected sequence (at the 3D structure resolution!)
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    beads_per_chromosome = [math.ceil(length/HiC_resolution) for length in chromosome_length]
    print(f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}")
    
    pdb_coordinates_df = pdb_coordinates.df["ATOM"]
    pdb_coordinates_df_output = pd.DataFrame()

    genome_fasta = SeqIO.to_dict(SeqIO.parse(fasta_name, "fasta"))

    for i, name in zip(range(len(chromosome_length)), chromosome_name):

        # Select one chromosome
        pdb_coordinates_chrom_x = pdb_coordinates_df[pdb_coordinates_df["residue_number"]==i+1]
        pdb_coordinates_chrom_x.reset_index(inplace=True, drop=True)

        # Calculate Euclydian distances between each pair of beads
        Euclydian_distances_after = np.sqrt((np.diff(np.array(pdb_coordinates_chrom_x["x_coord"]), axis=0))**2
                                                +(np.diff(np.array(pdb_coordinates_chrom_x["y_coord"]), axis=0))**2
                                                +(np.diff(np.array(pdb_coordinates_chrom_x["z_coord"]), axis=0))**2)
        Euclydian_distances_after = np.append(Euclydian_distances_after, [0])

        # Get chromosome ATCG sequence
        chromosome_sequence = str(genome_fasta[name].seq)

        # Select extremities of inverted contigs, i.e. beads further than 3*mean(distances) from the next bead.
        pdb_coordinates_chrom_x["distance"] = Euclydian_distances_after
        flipping_limits = pdb_coordinates_chrom_x["distance"]>3*np.mean(pdb_coordinates_chrom_x["distance"])
        flipping_limits_index = list(pdb_coordinates_chrom_x[flipping_limits].index)

        # Correct the beads order in the PDB and the bases order in the FASTA
        if flipping_limits_index != []:
            for j in range(0, len(flipping_limits_index), 2):
                
                # Flip the atoms corresponding to the inverted contigs
                print("Flipping atoms n°"+str(flipping_limits_index[j]+1)+" to n°"+str(flipping_limits_index[j+1]))
                flipped_atoms = pdb_coordinates_chrom_x.loc[flipping_limits_index[j]+1:flipping_limits_index[j+1],][::-1]
                pdb_coordinates_chrom_x = pd.concat([pdb_coordinates_chrom_x.iloc[:flipping_limits_index[j]+1,:], flipped_atoms.iloc[:,], pdb_coordinates_chrom_x.iloc[flipping_limits_index[j+1]+1:,]], axis=0)
                
                # Flip the bases corresponding th the inverted contigs
                print("Flipping bases n°"+str(flipping_limits_index[j]+1)+" to n°"+str(flipping_limits_index[j+1]+1))
                flipped_contig = chromosome_sequence[(flipping_limits_index[j]+1)*HiC_resolution:(flipping_limits_index[j+1]+1)*HiC_resolution]
                flipped_contig = flipped_contig[::-1]
                chromosome_sequence = chromosome_sequence[:(flipping_limits_index[j]+1)*HiC_resolution]+flipped_contig+chromosome_sequence[(flipping_limits_index[j+1]+1)*HiC_resolution:]
                genome_fasta[name].seq = Seq(chromosome_sequence)
        
        pdb_coordinates_chrom_x.reset_index(inplace=True, drop=True)
        pdb_coordinates_df_output = pd.concat([pdb_coordinates_df_output, pdb_coordinates_chrom_x.iloc[:,:-1]], axis=0)

    # Write corrected beads list and corrected bases sequence to PDB file and FASTA file
    with open(fasta_name_out, "w") as handle:
        SeqIO.write(genome_fasta.values(), handle, "fasta")

    pdb_coordinates_df_output.reset_index(inplace=True, drop=True)
    pdb_coordinates_df_output["line_idx"] = pdb_coordinates_df_output.index
    pdb_coordinates.df["ATOM"] = pdb_coordinates_df_output
    pdb_coordinates.to_pdb(path=pdb_name_out, records=None, gz=False, append_newline=True)
    print(f"Wrote {pdb_name_out} and {fasta_name_out}")


if __name__ == "__main__":
    ARGS = get_cli_arguments()

    # Read Fasta file and get chromosome length
    CHROMOSOME_LENGTH = extract_chromosome_length(ARGS.fasta)
    CHROMOSOME_NAME = extract_chromosome_name(ARGS.fasta)

    # Assign chromosome number
    flip_inverted_contig(ARGS.pdb, CHROMOSOME_LENGTH, CHROMOSOME_NAME, ARGS.fasta, ARGS.resolution, ARGS.output_pdb, ARGS.output_fasta)


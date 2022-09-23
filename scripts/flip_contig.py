"""Detect flipped contig from a PDB file containing a 3D genome structure.

Fix -- reinvert the flipped part of the genome -- in the 3D structure and the sequence.

This script requires:
- a PDB file containing the 3D genome structure,
- a fasta file containing the genome sequence,
- an Hi-C resolution.
"""

import argparse
import math

from Bio import SeqIO
from biopandas.pdb import PandasPdb
from Bio.Seq import Seq
import numpy as np
import pandas as pd


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
        "--fasta",
        action="store",
        type=str,
        help="Fasta file containing the sequence of the genome",
        required=True,
    )
    parser.add_argument(
        "--resolution",
        action="store",
        type=int,
        help="HiC resolution",
        required=True,
    )
    parser.add_argument(
        "--output-pdb",
        action="store",
        type=str,
        help="Output PDB file containing the fixed 3D structure of the genome",
        required=True,
    )
    parser.add_argument(
        "--output-fasta",
        action="store",
        type=str,
        help="Output FASTA file containing the fixed sequence of the genome",
        required=True,
    )
    return parser.parse_args()


def extract_chromosome_name_length(fasta_name):
    """Extract chromosome name and length from a FASTA file.

    Parameters
    ----------
    fasta_name : str
        Name of Fasta file containing the sequence of the genome
    
    Returns
    -------
    tuple
        List of chromosome names
        List of chromosome lengthes
    """
    chromosome_name_lst = []
    chromosome_length_lst = []
    with open(fasta_name, "r") as fasta_file:
        print(f"Reading {fasta_name}")
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id
            length = len(record.seq)
            print(f"Found chromosome {name} with {length} bases")
            chromosome_name_lst.append(name)
            chromosome_length_lst.append(length)
    return chromosome_name_lst, chromosome_length_lst



def find_inverted_contigs(pdb_name_in, chromosome_length, chromosome_name, fasta_name, HiC_resolution):
    """Detect inverted contigs.

    It uses the eucledian distance between adjacent beads in the 3D structure of the genome.

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
    """
    pdb_coordinates = PandasPdb().read_pdb(pdb_name_in)
    print(f"Number of beads read from structure: {pdb_coordinates.df['ATOM'].shape[0]}")

    beads_per_chromosome = [math.ceil(length/HiC_resolution) for length in chromosome_length]
    print(f"Number of beads deduced from sequence and HiC resolution: {sum(beads_per_chromosome)}")
    
    coordinates = pdb_coordinates.df["ATOM"]

    inverted_contigs = {}

    for chrom_index in coordinates["residue_number"].unique():
        print(f"Finding inverted contigs for chromosome {chrom_index}")
        
        # Select beads of one chromosome
        chrom_coordinates = coordinates.query(f"residue_number == {chrom_index}").reset_index(drop=True)
                
        # Compute Euclidean distances between bead n+1 and bead n
        euclidean_distances = np.sqrt( (np.diff(np.array(chrom_coordinates["x_coord"]), axis=0))**2
                                            +(np.diff(np.array(chrom_coordinates["y_coord"]), axis=0))**2
                                            +(np.diff(np.array(chrom_coordinates["z_coord"]), axis=0))**2
                                            )
        euclidean_distances = np.append(euclidean_distances, [0])

        mean_distance = np.mean(euclidean_distances)
        print(f"Mean distance between beads: {mean_distance:.2f}")
        print(f"Standard deviation of distances between beads: {np.std(euclidean_distances):.2f}")
        
        # Select extremities of inverted contigs
        # i.e. beads with distance beyond a given threshold of 3*mean(distances)
        chrom_coordinates = chrom_coordinates.assign(distance = euclidean_distances)
        beads_selection = chrom_coordinates["distance"]>3*mean_distance
        inversion_limits = chrom_coordinates.loc[beads_selection , "atom_number"].values
        print(inversion_limits)
        if len(inversion_limits)%2 != 0:
            print("WARNING: odd number of inversion limits found")

        inverted_contigs[chrom_index] = inversion_limits


def flip_inverted_contigs(pdb_name_in, chromosome_length, chromosome_name, fasta_name, HiC_resolution, pdb_name_out, fasta_name_out):
    """Flip contigs inverted in the genome assemblyaccording to the Euclidian distances between beads in the 3D genome structure.

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

        # Select the beads of one chromosome
        pdb_coordinates_chrom_x = pdb_coordinates_df[pdb_coordinates_df["residue_number"]==i+1]
        pdb_coordinates_chrom_x.reset_index(inplace=True, drop=True)

        # Calculate Euclydian distances between each pair of beads
        Euclydian_distances_after = np.sqrt((np.diff(np.array(pdb_coordinates_chrom_x["x_coord"]), axis=0))**2
                                                +(np.diff(np.array(pdb_coordinates_chrom_x["y_coord"]), axis=0))**2
                                                +(np.diff(np.array(pdb_coordinates_chrom_x["z_coord"]), axis=0))**2)
        Euclydian_distances_after = np.append(Euclydian_distances_after, [0])

        # Print mean and standard deviation of distances
        print("mean distance between beads : "+str(np.mean(Euclydian_distances_after)))
        print("Standard deviation of distances between beads : "+str(np.std(Euclydian_distances_after)))

        # Get chromosome ATCG sequence
        chromosome_sequence = str(genome_fasta[name].seq)

        # Select extremities of inverted contigs, i.e. beads further than 3*mean(distances) from the next bead (arbitrary)
        pdb_coordinates_chrom_x = pdb_coordinates_chrom_x.assign(distance = Euclydian_distances_after)
        flipping_limits = pdb_coordinates_chrom_x["distance"]>3*np.mean(pdb_coordinates_chrom_x["distance"])
        flipping_limits_index = list(pdb_coordinates_chrom_x[flipping_limits].index)

        # TO DO kill script if odd number of flipping_limits_index

        # Print number of inverted contigs, half the number of extremities
        print("number of inverted contigs : "+len(flipping_limits_index)/2)
        
        # Correct the beads order in the PDB and the bases order in the FASTA
        if flipping_limits_index != []:
            for j in range(0, len(flipping_limits_index), 2):

                # The beads inside an inverted contig are between a first extremity bead (excluded) and the next extremity bead (included)
                # The left extremity is excluded and the right one included because we look at the distance between the n and n+1 bead.
                print("Detected inverted contig between bead n°"+str(flipping_limits_index[j]+1)+" and bead n°"+str(flipping_limits_index[j+1]))
                
                # Flip the beads corresponding to the inverted contigs
                print("Flipping beads n°"+str(flipping_limits_index[j]+1)+" to n°"+str(flipping_limits_index[j+1]))
                flipped_beads = pdb_coordinates_chrom_x.loc[flipping_limits_index[j]+1:flipping_limits_index[j+1],][::-1]
                pdb_coordinates_chrom_x = pd.concat([pdb_coordinates_chrom_x.iloc[:flipping_limits_index[j]+1,:], flipped_beads.iloc[:,], pdb_coordinates_chrom_x.iloc[flipping_limits_index[j+1]+1:,]], axis=0)
                
                # Flip the bases corresponding th the inverted contigs
                print("Flipping bases n°"+str((flipping_limits_index[j]+1)*HiC_resolution)+" to n°"+str((flipping_limits_index[j+1]+1)*HiC_resolution))
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
    print(ARGS)
    # Read Fasta file and extract chromosome name and length
    CHROMOSOME_NAME, CHROMOSOME_LENGTH = extract_chromosome_name_length(ARGS.fasta)

    # Assign chromosome number
    find_inverted_contigs(ARGS.pdb, CHROMOSOME_LENGTH, CHROMOSOME_NAME, ARGS.fasta, ARGS.resolution)
    flip_inverted_contigs(ARGS.pdb, CHROMOSOME_LENGTH, CHROMOSOME_NAME, ARGS.fasta, ARGS.resolution, ARGS.output_pdb, ARGS.output_fasta)


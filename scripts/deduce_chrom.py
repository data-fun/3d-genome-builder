from pysam import FastaFile
from biopandas.pdb import PandasPdb

pdb_name_1 = "./sample-B_50000_poster.pdb"
fasta_name = "./Neurospora_nc14.fasta"
HiC_resolution = 50000

def deduce_chrom(pdb_name, fasta_name, HiC_resolution, output_file):
    """Assign each point to a chromosome, depending on the resolution of the model and the size of the chromosomes."""
    
    # get chromosomes sizes from fasta file
    genome_sequence = FastaFile(fasta_name)
    chrom_lengths = genome_sequence.lengths
    
    # read .pdb file
    atoms_coordinates = PandasPdb().read_pdb(pdb_name)
    
    # count atoms associated with each chromosome according to the resolution
    nb_atoms = []
    for i in range(len(chrom_lengths)):
        nb_atoms_chrom_x = len(range(0, chrom_lengths[i], HiC_resolution))
        atoms_coordinates.df["ATOM"].loc[sum(nb_atoms):sum(nb_atoms)+nb_atoms_chrom_x-1,"residue_number"] = str(i+1)
        nb_atoms.append(nb_atoms_chrom_x)
    
    atoms_coordinates.df["ATOM"]["residue_name"] = "CHR"

    # return .pdb with chromosomes numbers in the column 5 nd 7 : "CHR" (or "PLA" ?) with the number
    atoms_coordinates.to_pdb(path=output_file, records=None, gz=False, append_newline=True)



deduce_chrom(pdb_name_1, fasta_name, HiC_resolution, "./sampleB_50000_formalized.pdb")
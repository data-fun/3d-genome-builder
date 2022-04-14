
"""
main_config.txt ?
species_name = 
HiC fastq SRA_IDs = []
restriction_enzyme_name = [] #one or more
"""

#CONDA ENV 3D_models

#loop on nb of SRA_ID ? Put (compressed ?) files into separated folders (./species_name/raw_data/sample_X)
fastq-dump SRR16761089 --split-files
fastq-dump SRR14362684 --split-files
fastq-dump SRR16761090 --split-files
fastq-dump SRR16761091 --split-files
fastq-dump SRR16761092 --split-files

#merge_fastq -fp1 SRR16761092_1.fastq -fp1 SRR16761089_1.fastq -fp1 SRR14362684_1.fastq -fp1 SRR16761090_1.fastq -fp1 SRR16761091_1.fastq -fp2 SRR16761092_2.fastq -fp2 SRR16761089_2.fastq -fp2 SRR14362684_2.fastq -fp2 SRR16761090_2.fastq -fp2 SRR16761091_2.fastq

#Get genome .fasta file
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173593/suppl/GSE173593%5Fnc14%5Fgenome%2Dseq%2Dname%2Efasta%2Egz
#rename chromosome, delete scafold

#Get info on chromosomes sizes from the fasta file ? + delete scafolds size --> a .txt file with two columns, chrom_name and size.
#Put the .txt file "species_name_chrom_sizes.txt" in ./species_name/HiC-Pro_config_files

#Generate the list of restriction fragments after genome digestion
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r dpnii T^TAA -o ./Neurospora/HiC-Pro_config_files/Neurospora_DpnII-MseI ./Neurospora/Neurospora_genome.fasta

#Generate Bowtie2 index
bowtie2-build ./Neurospora_genome.fasta neurospora
#Put result files in species_name/HiC-Pro_config_files/index

#Run HiC-Pro
HiC-Pro -i ./raw_data -o ./output -c ./HiC-Pro_config_files/config-hicpro.txt

#Merge samples 1 to 5 by reruning HiC-Pro on .validpairs in merge_data folder https://github.com/nservant/HiC-Pro/issues/121
HiC-Pro -i ./merged_samples -o ./output -c ./HiC-Pro_config_files/config-hicpro.txt -s merge_persample -s build_contact_maps -s ice_norm


#download HiCPlotter
https://gitee.com/simonjyoung/HiCPlotter/blob/master/HiCPlotter.py

#Run HiCPlotter
#10K resolution
python ./HiCPlotter/HiCPlotter.py -f output/hic_results/matrix/merged_samples/iced/10000/merged_samples_10000_iced.matrix -o species_name -r 10000 -tri 1 -bed output/hic_results/matrix/sample_merge/raw/10000/sample_merge_10000_abs.bed -n TEST -wg 1 -chr NC_026507.1 -hmc 4
#50K resolution
python ./HiCPlotter/HiCPlotter.py -f output/hic_results/matrix/merged_samples/iced/50000/merged_samples_50000_iced.matrix -o species_name -r 50000 -tri 1 -bed output/hic_results/matrix/sample_merge/raw/50000/sample_merge_50000_abs.bed -n TEST -wg 1 -chr NC_026507.1 -hmc 4

#####################

## Convert to dense format
/home/thibault/HiC-Pro-master/bin/utils/sparseToDense.py -b output/hic_results/matrix/sample_merge/raw/40000/sample_merge_40000_abs.bed output/hic_results/matrix/sample_merge/iced/40000/sample_merge_40000_iced.matrix
#put result into dense_matrix folder

#convert to numpy array
"""PYTHON"""
import numpy as np
raw_matrix_50K = np.loadtxt("./dense_matrix/merge_samples_50000_dense.matrix")
np.save("counts_matrix_50K.npy", raw_matrix_50K)
"""PYTHON"""
#Put result into ./pastis/data
#create a species_name_structure file with only chrom sizes (one by rows) and put it in pastis/files
pastis-pm2 <work_folder>

##create csv from pdb : rename .pdb file to .txt and .pdb.pdb file to .pdb then convert this one to .csv (to get tab separated values and not multiple spaces separated values)
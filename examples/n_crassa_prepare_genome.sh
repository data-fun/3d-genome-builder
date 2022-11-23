# Create target folder
target_folder="3DGB_n_crassa_WT"
mkdir -p "${3DGB_n_crassa}"

# Download genome sequence
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173593/suppl/GSE173593%5Fnc14%5Fgenome%2Dseq%2Dname%2Efasta%2Egz -O "${target_folder}/genome.fasta.gz"
gunzip "${target_folder}/genome.fasta.gz"

# Rename chromosome
sed -i "s/Supercontig_14./chr/" "${target_folder}/genome.fasta"

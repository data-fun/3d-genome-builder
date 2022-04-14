from Bio import SeqIO
with open(snakemake.input[0], "r") as fasta_file, open(snakemake.output[0], "w") as size_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"{record.id}: {len(record.seq)} bases")
        size_file.write(f"{record.id}\t{len(record.seq)}\n")
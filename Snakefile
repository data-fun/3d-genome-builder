print(config)

WORKDIR = f"3d_{config['organism_name_short'].lower()}"

rule all:
    input:
        #expand("{workdir}/raw_data/{sra_id}/{sra_id}_{paired}.fastq.gz", 
        #       workdir=WORKDIR, sra_id=config['sra_ids'], paired=[1])
        expand("{workdir}/HiC-Pro/chromosome_sizes.txt", workdir=WORKDIR),
        expand("{workdir}/HiC-Pro/restriction_sites.txt", workdir=WORKDIR),
        expand("{workdir}/HiC-Pro/bowtie2_index/genome.1.bt2", workdir=WORKDIR)


# fasterq-dump documentation:
# https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
rule download_fastq:
    output:
        "{WORKDIR}/raw_data/{sra_id}/{sra_id}_1.fastq"
    params:
        "--progress"
    threads:
        4
    log:
        "{WORKDIR}/logs/fasterq-dump_{sra_id}.log"
    message:
        "Downloading {wildcards.sra_id} files"
    shell: 
        "mkdir -p {WORKDIR}/raw_data/{wildcards.sra_id} && "
        "fasterq-dump {wildcards.sra_id} --threads {threads} {params} --outdir {WORKDIR}/raw_data/{wildcards.sra_id}/ "


rule compress_fastq:
    input:
        "{WORKDIR}/raw_data/{sra_id}/{sra_id}_1.fastq"
    output:
        "{WORKDIR}/raw_data/{sra_id}/{sra_id}_1.fastq.gz"
    threads:
        4
    log:
        "{WORKDIR}/logs/gzip_{sra_id}.log"
    message:
        "Compressing {input}"
    shell: 
        # pigz is the parallel implementation of gzip
        "pigz -p {threads} {WORKDIR}/raw_data/{wildcards.sra_id}/*"


rule calculate_chromosome_sizes:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/chromosome_sizes.txt"
    message:
        "Calculating chromosome sizes"
    script:
        "scripts/calculate_chromosome_sizes.py"


    
rule hicpro:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/restriction_sites.txt"
    container:
        "images/hicpro.img"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py "
        "-r {config[hicpro_restriction_sites]} "
        "-o {output} "
        "{input}"


rule build_genome_index:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/bowtie2_index/genome.1.bt2"
    container:
        "images/hicpro.img"
    shell:
        "bowtie2-build {input} {WORKDIR}/HiC-Pro/bowtie2_index/genome"

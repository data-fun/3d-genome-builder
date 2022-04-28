configfile: "config.yml"
WORKDIR = f"3d_{config['organism_name_short'].lower()}"
workdir: f"{WORKDIR}"

rule all:
    input:
        expand("fastq_files/{sra_id}/{sra_id}_R{paired}.fastq.gz",
            sra_id=config['sra_ids'], paired=[1,2]),
        "HiC-Pro/config.txt",
        expand("HiC-Pro/output/logs/{sra_id}/build_raw_maps.log",
            sra_id=config['sra_ids'])



# fasterq-dump documentation:
# https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
rule download_fastq_files:
    output:
        R1="fastq_files/{sra_id}/{sra_id}_1.fastq",
        R2="fastq_files/{sra_id}/{sra_id}_2.fastq"
    params:
        outdir="fastq_files/{sra_id}/",
        progress="--progress",
        loglevel="info"
    threads:
        2
    log:
        "logs/fasterq-dump_{sra_id}.log"
    message:
        "Downloading {wildcards.sra_id} files"
    shell: 
        "fasterq-dump {wildcards.sra_id} --threads {threads} --log-level {params.loglevel} "
        "{params.progress} --outdir {params.outdir} 2>&1 >{log}"


rule rename_paired_fastq_files:
    input:
        "fastq_files/{sra_id}/{sra_id}_{paired}.fastq"
    output:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq"
    message:
        "Renamming {wildcards.sra_id} files"
    shell: 
        "mv {input} {output}"


rule compress_fastq_files:
    input:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq"
    output:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq.gz"
    threads:
        4
    log:
        "logs/pigz_{sra_id}_R{paired}.log"
    message:
        "Compressing {input}"
    shell: 
        # pigz is the parallel implementation of gzip
        "pigz -p {threads} {input}"


rule calculate_chromosome_sizes:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/chromosome_sizes.txt"
    message:
        "Calculating chromosome sizes"
    script:
        "scripts/calculate_chromosome_sizes.py"


rule digest_genome_HiC_Pro:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/restriction_sites.txt"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py "
        "-r {config[hicpro_restriction_sites]} "
        "-o {output} "
        "{input}"


rule build_genome_index_HiC_Pro:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/bowtie2_index/genome.1.bt2"
    log:
        "logs/build_genome_index.log"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    shell:
        "bowtie2-build {input} HiC-Pro/bowtie2_index/genome 2>&1 | tee {log}"


rule create_HiC_Pro_config:
    input:
        template="../templates/HiC-Pro_config.template",
        chromosome_sizes="HiC-Pro/chromosome_sizes.txt",
        genome_fragment="HiC-Pro/restriction_sites.txt"
    output:
        "HiC-Pro/config.txt"
    params:
        genome_index_path="HiC-Pro/bowtie2_index", 
    script:
        "scripts/create_HiC_Pro_config.py"


rule run_HiC_Pro:
    input:
        config="HiC-Pro/config.txt",
        genome_index="HiC-Pro/bowtie2_index/genome.1.bt2"
    output:
        expand("HiC-Pro/output/logs/{sra_id}/build_raw_maps.log",
            sra_id=config['sra_ids'])
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    threads: 8
    shell:
        "HiC-Pro -i fastq_files -o HiC-Pro/output -c {input.config}"

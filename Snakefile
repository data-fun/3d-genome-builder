#print(config)

WORKDIR = f"3d_{config['organism_name_short'].lower()}"

rule all:
    input:
        expand("{workdir}/fastq_files/{sra_id}/{sra_id}_R{paired}.fastq.gz",
               workdir=WORKDIR, sra_id=config['sra_ids'], paired=[1]),
        expand("{workdir}/HiC-Pro/chromosome_sizes.txt", workdir=WORKDIR),
        expand("{workdir}/HiC-Pro/restriction_sites.txt", workdir=WORKDIR),
        expand("{workdir}/HiC-Pro/bowtie2_index/genome.1.bt2", workdir=WORKDIR),
        expand("{workdir}/HiC-Pro/config.txt", workdir=WORKDIR)


# fasterq-dump documentation:
# https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
rule download_fastq_files:
    output:
        R1="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_1.fastq",
        R2="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_2.fastq"
    params:
        outdir="{WORKDIR}/fastq_files/{sra_id}/",
        progress="--progress",
        loglevel="info"
    threads:
        2
    log:
        "{WORKDIR}/logs/fasterq-dump_{sra_id}.log"
    message:
        "Downloading {wildcards.sra_id} files"
    shell: 
        "fasterq-dump {wildcards.sra_id} --threads {threads} --log-level {params.loglevel} "
        "{params.progress} --outdir {params.outdir} 2>&1 >{log}"


rule rename_paired_fastq_files:
    input:
        R1="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_1.fastq",
        R2="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_2.fastq"
    output:
        R1="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R1.fastq",
        R2="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R2.fastq"
    message:
        "Renamming {wildcards.sra_id} files"
    shell: 
        "mv {input.R1} {output.R1} &&"
        "mv {input.R2} {output.R2} "


rule compress_fastq_files:
    input:
        R1="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R1.fastq",
        R2="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R2.fastq",
    output:
        R1="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R1.fastq.gz",
        R2="{WORKDIR}/fastq_files/{sra_id}/{sra_id}_R2.fastq.gz"
    threads:
        4
    log:
        "{WORKDIR}/logs/pigz_{sra_id}.log"
    message:
        "Compressing {input}"
    shell: 
        # pigz is the parallel implementation of gzip
        "pigz -p {threads} {input}"


rule calculate_chromosome_sizes:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/chromosome_sizes.txt"
    message:
        "Calculating chromosome sizes"
    script:
        "scripts/calculate_chromosome_sizes.py"


    
rule digest_genome_HiC_Pro:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/restriction_sites.txt"
    container:
        "images/hicpro.sif"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py "
        "-r {config[hicpro_restriction_sites]} "
        "-o {output} "
        "{input}"


rule build_genome_index_HiC_Pro:
    input:
        "{WORKDIR}/genome.fasta"
    output:
        "{WORKDIR}/HiC-Pro/bowtie2_index/genome.1.bt2"
    log:
        "{WORKDIR}/logs/build_genome_index.log"
    container:
        "images/hicpro.sif"
    shell:
        "bowtie2-build {input} {WORKDIR}/HiC-Pro/bowtie2_index/genome 2>&1 | tee {log}"


rule create_HiC_Pro_config:
    input:
        template="templates/HiC-Pro_config.template",
        chromosome_sizes="{WORKDIR}/HiC-Pro/chromosome_sizes.txt",
        genome_fragment="{WORKDIR}/HiC-Pro/restriction_sites.txt"
    output:
        "{WORKDIR}/HiC-Pro/config.txt"
    params:
        genome_index_path="{WORKDIR}/HiC-Pro/bowtie2_index", 
    script:
        "scripts/create_HiC_Pro_config.py"


rule run_HiC_Pro:
    input:
        config="3d_n_crassa/HiC-Pro/config.txt"
    container:
        "images/hicpro.sif"
    shell:
        "HiC-Pro -i ./3d_n_crassa/raw_data -o ./3d_n_crassa/output -c {input.config}"
configfile: "config.yml"

ORGANISM_NAME = f"{config['organism_name_short'].lower()}"
WORKDIR = f"3d_{ORGANISM_NAME}"
workdir: f"{WORKDIR}"

rule all:
    input:
        expand("HiCPlotter/hicplotter_{resolution}.ok", resolution=config["hicpro_resolutions"])



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
    conda:
        "envs/workflow.yml"
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
    conda:
        "envs/workflow.yml"
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
    conda:
        "envs/workflow.yml"
    script:
        "scripts/calculate_chromosome_sizes.py"


# Generate the list of restriction fragments
# with a tool provided with HiC-Pro
rule digest_genome:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/restriction_sites.txt"
    message:
        "Generating list of restriction fragments"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py "
        "-r {config[hicpro_restriction_sites]} "
        "-o {output} "
        "{input}"


# Build genome index
# with bowtie provided with HiC-Pro
rule build_genome_index:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/bowtie2_index/genome.1.bt2"
    log:
        "logs/build_genome_index.log"
    message:
        "Building genome index"
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
        resolutions=config["hicpro_resolutions"]
    message:
        "Building HiC-Pro configuration file"
    conda:
        "envs/workflow.yml"
    script:
        "scripts/create_HiC_Pro_config.py"


rule run_HiC_Pro:
    input:
        config="HiC-Pro/config.txt",
        genome_index="HiC-Pro/bowtie2_index/genome.1.bt2",
        fastq_files=expand("fastq_files/{sra_id}/{sra_id}_R{paired}.fastq.gz",
                    sra_id=config["sra_ids"], paired=[1,2])
    output:
        expand("HiC-Pro/output/hic_results/data/{sra_id}/{sra_id}_genome.bwt2pairs.validPairs",
        sra_id=config["sra_ids"])
    message:
        "Running HiC-Pro"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    threads: 8
    shell:
        "HiC-Pro -i fastq_files -o HiC-Pro/output -c {input.config}"


# Merge validPairs produced by HiC-Prot
# See https://github.com/nservant/HiC-Pro/issues/121
rule copy_HiC_Pro_valid_pairs:
    input:
        "HiC-Pro/output/hic_results/data/{sra_id}/{sra_id}_genome.bwt2pairs.validPairs"
    output:
        "HiC-Pro/merged_samples/merge/{sra_id}_genome.bwt2pairs.validPairs"
    message:
        "Merging valid pairs"
    shell:
        "cp {input} {output}"


rule run_HiC_Pro_on_valid_pairs:
    input:
        config="HiC-Pro/config.txt",
        valid_pairs=expand("HiC-Pro/merged_samples/merge/{sra_id}_genome.bwt2pairs.validPairs",
                    sra_id=config['sra_ids'])
    output:
        matrix=expand("HiC-Pro/merged_output/hic_results/matrix/merge/iced/{resolution}/merge_{resolution}_iced.matrix",
               resolution=config["hicpro_resolutions"]),
        bed=expand("HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed",
            resolution=config["hicpro_resolutions"])
    message:
        "Running HiC-Pro on valid pairs"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    threads: 8
    shell:
        "HiC-Pro -i HiC-Pro/merged_samples -o HiC-Pro/merged_output -c {input.config} -s merge_persample -s build_contact_maps -s ice_norm"


rule run_HiCPlotter:
    input:
        matrix="HiC-Pro/merged_output/hic_results/matrix/merge/iced/{resolution}/merge_{resolution}_iced.matrix",
        bed="HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed"
    output:
        "HiCPlotter/hicplotter_{resolution}.ok"
    message:
        "Running HiCPlotter"
    conda:
        "envs/workflow.yml"
    shell:
        "python ../scripts/HiCPlotter.py "
        "-f {input.matrix} "
        "-bed {input.bed} "
        f"-o HiCPlotter/{ORGANISM_NAME} "
        "-r {wildcards.resolution} "
        "-tri 1  "         # input file is from HiC-Pro pipeline
        "-n TEST "
        "-wg 1 "           # plotting whole genome interactions
        "-chr chr1 "       # chromosome to be plotted
        "-hmc 4 && "
        "touch {output}"   # dummy output files

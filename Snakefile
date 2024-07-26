import os

ORGANISM_NAME = f"{config['organism'].replace(' ', '_')}"
WORKING_DIR = f"{config['workdir'].replace(' ', '_')}"

# Verify the present of 'hicpro_ligation_site'
# and set a default value if not.
if "hicpro_ligation_site" not in config:
    config["hicpro_ligation_site"] = "no_ligation_site"


workdir: WORKING_DIR

#  Main rule to generate all result files.
rule all:
    input:
        "logs/environment_conda.log",
        "logs/environment_singularity.log",
        "HiC-Pro/config.txt",
        expand("contact_maps/contact_map_{resolution}.png", resolution=config["hicpro_resolutions"]),
        expand("pastis/structure_{resolution}.pdb", resolution=config["pastis_resolutions"]),
        expand("structure/{resolution}/structure_with_chr.pdb", resolution=config["pastis_resolutions"]),
        expand("structure/{resolution}/structure_completed.pdb", resolution=config["pastis_resolutions"]),
        expand("structure/{resolution}/structure_cleaned.pdb", resolution=config["pastis_resolutions"]),
        expand("structure/{resolution}/structure_cleaned.g3d", resolution=config["pastis_resolutions"]),


# Extract versions of Python modules and tools from conda environment.
rule describe_conda_env:
    priority: 50
    output:
        "logs/environment_conda.log"
    message:
        "Extracting information from conda environment"
    conda:
        "envs/workflow.yml"
    shell:
        "python ../scripts/describe_conda_env.py >{output} &&"
        "fasterq-dump --version >>{output} &&"
        "pigz --version >>{output} "


# Extract versions of tools from singularity container.
rule describe_singularity_env:
    priority: 51
    output:
        "logs/environment_singularity.log"
    message:
        "Extracting information from singularity container"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py --help >{output} && "
        "printf '\n\n\n' >>{output} && "
        "bowtie2-build --version >>{output} &&"
        "printf '\n\n\n' >>{output} && "
        "LC_ALL=C; HiC-Pro --version >>{output} "



# Donwload fastq files.
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
    message:
        "Downloading {wildcards.sra_id} files"
    conda:
        "envs/workflow.yml"
    log:
        "logs/fasterq-dump_{sra_id}.log"
    shell: 
        "fasterq-dump {wildcards.sra_id} --threads {threads} --log-level {params.loglevel} "
        "{params.progress} --outdir {params.outdir} >{log} 2>&1 "


# Rename fastq files to be compatible with HiC-Pro.
rule rename_paired_fastq_files:
    input:
        "fastq_files/{sra_id}/{sra_id}_{paired}.fastq"
    output:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq"
    message:
        "Renamming {wildcards.sra_id} files"
    shell: 
        "mv {input} {output}"


# Compress fastq files to gain some space.
rule compress_fastq_files:
    input:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq"
    output:
        "fastq_files/{sra_id}/{sra_id}_R{paired}.fastq.gz"
    threads:
        4
    message:
        "Compressing {input}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/pigz_{sra_id}_R{paired}.log"
    shell: 
        # pigz is the parallel implementation of gzip
        "pigz -p {threads} {input} >{log} 2>&1"


# Calculate chromosome sizes from the genome sequence.
rule calculate_chromosome_sizes:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/chromosome_sizes.txt"
    message:
        "Calculating chromosome sizes"
    conda:
        "envs/workflow.yml"
    log:
        "logs/calculate_chromosome_sizes.log"
    shell:
        "python ../scripts/calculate_chromosome_sizes.py "
        "--fasta {input} "
        "--output {output} "
        ">{log} 2>&1 "


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
    log:
        "logs/digest_genome.log"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py "
        "-r {config[hicpro_restriction_sites]} "
        "-o {output} "
        "{input} "
        ">{log} 2>&1 "


# Build genome index
# with bowtie provided with HiC-Pro.
rule build_genome_index:
    input:
        "genome.fasta"
    output:
        "HiC-Pro/bowtie2_index/genome.1.bt2"
    message:
        "Building genome index"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    log:
        "logs/build_genome_index.log"
    shell:
        "bowtie2-build {input} HiC-Pro/bowtie2_index/genome >{log} 2>&1"


# Create HiC-Pro configuration file.
# See https://nservant.github.io/HiC-Pro/MANUAL.html#setting-the-configuration-file
rule create_HiC_Pro_config:
    input:
        template="../templates/HiC-Pro_config.template",
        chromosome_sizes="HiC-Pro/chromosome_sizes.txt",
        genome_fragment="HiC-Pro/restriction_sites.txt"
    output:
        "HiC-Pro/config.txt"
    params:
        genome_index_path="HiC-Pro/bowtie2_index",
        ligation_site=config["hicpro_ligation_site"],
        resolutions=config["hicpro_resolutions"]
    message:
        "Building HiC-Pro configuration file"
    conda:
        "envs/workflow.yml"
    log:
        "logs/create_HiC_Pro_config.log"
    shell:
        "python ../scripts/create_HiC_Pro_config.py "
        "--template {input.template} "
        "--chromosome-sizes {input.chromosome_sizes} "
        "--genome-fragment {input.genome_fragment} "
        "--ligation-site {params.ligation_site} "
        "--genome-index-path {params.genome_index_path} "
        "--resolutions {params.resolutions} "
        "--output {output} "
        ">{log} 2>&1"


# Run HiC-Pro
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
    log:
        "logs/run_HiC_Pro.log"
    shell:
        "LC_ALL=C; echo 'y' | HiC-Pro "
        "-i fastq_files "
        "-o HiC-Pro/output "
        "-c {input.config} "
        ">{log} 2>&1 "
    # LC_ALL=C prevents warning messages on locale settings
    # /usr/bin/bash: warning: setlocale: LC_ALL: cannot change locale (fr_FR.UTF-8)
    # echo 'y' automatically valids answer to
    # [...]/HiC-Pro/output folder alreads exists. Do you want to overwrite it ? (y/n) [n] :


# Merge validPairs produced by HiC-Pro
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


# Run HiC-Pro on merged valid pairs.
rule run_HiC_Pro_on_valid_pairs:
    input:
        config="HiC-Pro/config.txt",
        valid_pairs=expand("HiC-Pro/merged_samples/merge/{sra_id}_genome.bwt2pairs.validPairs",
                           sra_id=config["sra_ids"])
    output:
        matrix_iced=expand("HiC-Pro/merged_output/hic_results/matrix/merge/iced/{resolution}/merge_{resolution}_iced.matrix",
                      resolution=config["hicpro_resolutions"]),
        matrix_raw=expand("HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}.matrix",
                   resolution=config["hicpro_resolutions"]),
        bed=expand("HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed",
                   resolution=config["hicpro_resolutions"])
    message:
        "Running HiC-Pro on valid pairs"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    threads: 8
    log:
        "logs/run_HiC_Pro_on_valid_pairs.log"
    shell:
        "HiC-Pro "
        "-i HiC-Pro/merged_samples "
        "-o HiC-Pro/merged_output "
        "-c {input.config} "
        "-s merge_persample "
        "-s build_contact_maps "
        "-s ice_norm "
        ">{log} 2>&1 "


rule convert_iced_matrix_sparse_to_dense:
    input:
        bed="HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed",
        matrix="HiC-Pro/merged_output/hic_results/matrix/merge/iced/{resolution}/merge_{resolution}_iced.matrix"
    output:
        "dense_matrix/iced/merge_{resolution}_dense.matrix"
    message:
        "Converting sparse to dense matrix (iced)"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    log:
        "logs/convert_iced_matrix_{resolution}.log"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/sparseToDense.py "
        "--bins {input.bed} "
        "--output {output} "
        "{input.matrix} "
        ">{log} 2>&1 "


rule convert_raw_matrix_sparse_to_dense:
    input:
        bed="HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed",
        matrix="HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}.matrix"
    output:
        "dense_matrix/raw/merge_{resolution}_dense.matrix"
    message:
        "Converting sparse to dense matrix (raw)"
    container:
        "../images/hicpro_3.1.0_ubuntu.img"
    log:
        "logs/convert_raw_matrix_{resolution}.log"
    shell:
        "/usr/local/bin/HiC-Pro_3.1.0/bin/utils/sparseToDense.py "
        "--bins {input.bed} "
        "--output {output} "
        "{input.matrix} "
        ">{log} 2>&1 "


rule build_contact_map:
    input:
        "dense_matrix/iced/merge_{resolution}_dense.matrix"
    output:
        "contact_maps/contact_map_{resolution}.png"
    message:
        "Assembling contact maps"
    conda:
        "envs/workflow.yml"
    log:
        "logs/build_contact_map_{resolution}.log"
    shell:
        "python ../scripts/build_contact_maps.py "
        "--contacts {input} "
        "--map {output} "
        ">{log} 2>&1 "


rule run_pastis_nb:
    input:
        matrix="dense_matrix/raw/merge_{resolution}_dense.matrix",
        bed="HiC-Pro/merged_output/hic_results/matrix/merge/raw/{resolution}/merge_{resolution}_abs.bed"
    output:
        "pastis/structure_{resolution}.pdb"
    message:
        "Running Pastis NB at resolution {wildcards.resolution}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/run_pastis_nb_{resolution}.log"
    shell:
        "python ../scripts/infer_structures_nb.py "
        "--matrix {input.matrix} "
        "--bed {input.bed} "
        "--output {output} "
        ">{log} 2>&1 "


rule assign_chromosomes:
    input:
        structure="pastis/structure_{resolution}.pdb",
        sequence="genome.fasta"
    output:
        "structure/{resolution}/structure_with_chr.pdb"
    message:
        "Assigning chromosomes to the 3D structure at resolution {wildcards.resolution}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/assign_chromosomes_{resolution}.log"
    shell:
        "python ../scripts/assign_chromosomes.py "
        "--pdb {input.structure} "
        "--fasta {input.sequence} "
        "--resolution {wildcards.resolution} "
        "--output {output} "
        ">{log} 2>&1 "
        

rule verify_inverted_contigs:
    input:
        structure="structure/{resolution}/structure_with_chr.pdb",
        sequence="genome.fasta"
    output:
        structure="structure/{resolution}/structure_verified_contigs.pdb",
        sequence="sequence/{resolution}/genome_verified_contigs.fasta"
    params:
        verify=str(config.get("verify_contigs", False))
    message:
        "Fix inverted contigs (if needed) in the 3D structure at resolution {wildcards.resolution}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/verify_inverted_contigs_{resolution}.log"
    shell:
        "python ../scripts/verify_inverted_contigs.py "
        "--run {params.verify} "
        "--pdb {input.structure} "
        "--fasta {input.sequence} "
        "--resolution {wildcards.resolution} "
        "--output-pdb {output.structure} "
        "--output-fasta {output.sequence} "
        ">{log} 2>&1 "


rule add_missing_beads:
    input:
        "structure/{resolution}/structure_verified_contigs.pdb"
    output:
        "structure/{resolution}/structure_completed.pdb"
    message:
        "Adding missing beads in the 3D structure at resolution {wildcards.resolution}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/add_missing_beads_{resolution}.log"
    shell:
        "python ../scripts/add_missing_beads.py "
        "--input-pdb {input} "
        "--output-pdb {output} "
        ">{log} 2>&1 "


rule delete_outlier_beads:
    input:
        "structure/{resolution}/structure_completed.pdb"
    output:
        "structure/{resolution}/structure_cleaned.pdb"
    message:
        "Deleting outliers beads in the 3D structure at resolution {wildcards.resolution}"
    conda:
        "envs/workflow.yml"
    log:
        "logs/delete_outlier_beads_{resolution}.log"
    shell:
        "python ../scripts/delete_outlier_beads.py "
        "--input-pdb {input} "
        "--output-pdb {output} "
        ">{log} 2>&1 "


rule convert_to_g3d:
    input:
        structure="structure/{resolution}/structure_cleaned.pdb",
        sequence="genome.fasta"
    output:
        "structure/{resolution}/structure_cleaned.g3d"
    message:
        "Converting to g3d format"
    conda:
        "envs/workflow.yml"
    log:
        "logs/convert_to_g3d_{resolution}.log"
    shell:
        "python ../scripts/convert_to_g3d.py "
        "--pdb {input.structure} "
        "--fasta {input.sequence} "
        "--resolution {wildcards.resolution} "
        "--output {output} "
        ">{log} 2>&1 "


onsuccess:
    print("WORKFLOW COMPLETED SUCCESSFULLY!")
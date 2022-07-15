workdir: f"{config['workdir'].replace(' ', '_')}"


rule all:
    input:
        expand("structure/{resolution}/structure_with_chr.pdb", resolution=config["pastis_resolutions"]),
        expand("structure/{resolution}/structure_completed.pdb", resolution=config["pastis_resolutions"]),


rule assign_chromosomes:
    input:
        structure="pastis/{resolution}/structure.pdb",
        sequence="genome.fasta"
    output:
        "structure/{resolution}/structure_with_chr.pdb"
    message:
        "Assigning chromosomes"
    conda:
        "envs/workflow_structure.yml"
    shell:
        "python ../scripts/assign_chromosomes.py "
        "--pdb {input.structure} "
        "--fasta {input.sequence} "
        "--resolution {wildcards.resolution} "
        "--output {output}"
        

rule interpolate_missing_coordinates:
    input:
        "structure/{resolution}/structure_with_chr.pdb"
    output:
        "structure/{resolution}/structure_completed.pdb"
    message:
        "Fixing missing coordinates"
    conda:
        "envs/workflow_structure.yml"
    shell:
        "python ../scripts/interpolate_missing_coordinates.py "
        "--pdb {input} "
        "--output {output}"
# 3D genome builder

## Download this repository

```bash
git clone https://github.com/data-fun/3d-genome-builder.git
cd 3d-genome-builder
```

## Install dependencies

### Singularity

Download last version [here](https://github.com/apptainer/singularity/releases)

Install package:

```bash
sudo apt install -y ./singularity-container_3.8.7_amd64.deb

Verify version:

```
$ singularity --version
singularity version 3.8.7
```

### HiCPlotter

```bash
wget https://gitee.com/simonjyoung/HiCPlotter/raw/master/HiCPlotter.py -P scripts
```

### Conda environment

Install conda.

Install mamba:

```bash
conda install mamba -n base -c conda-forge
```

Create conda environment and install dependendies:

```bash
mamba env create -f binder/environment.yml
```

Load conda environment:

```bash
conda activate 3d-genome-builder
```

## Download  HiC-Pro Singularity image


```bash
wget --ciphers=DEFAULT:@SECLEVEL=1 https://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_3.1.0_ubuntu.img -P images
```

Verify HiC-Pro version with:

```bash
$ singularity exec images/hicpro.sif HiC-Pro --version
[...]
HiC-Pro version 3.1.0
```

and bowtie2 version:

```bash
$ singularity exec images/hicpro.sif bowtie2 --version  2>/dev/null | head -n 1
/usr/local/conda/envs/hicpro/bin/bowtie2-align-s version 2.4.4
```


## Add parameters

Create and edit the configuration file in yaml format. See for instance `config.yml`

## Add reference genome (mandatory)

The reference genome fasta file must be located in `3d_genome_n_crassa/genome.fasta` where `3d_genome_n_crassa` is the name of the working directory as specified in the config file `config.yml`.

## Add FASTQ files (optional)

If you already have fastq files stored locally or fastq files are not available on GEO or SRA, you can use these files providing they are in the proper directory structure:

```
3d_genome_n_crassa/
├── fastq_files
│   ├── SRR2105869
│   │   ├── SRR2105869_R1.fastq.gz
│   │   └── SRR2105869_R2.fastq.gz
│   ├── SRR2105870
│   │   ├── SRR2105870_R1.fastq.gz
│   │   └── SRR2105870_R2.fastq.gz
│   ├── SRR2105871
│   │   ├── SRR2105871_R1.fastq.gz
│   │   └── SRR2105871_R2.fastq.gz
│   └── SRR2105872
│       ├── SRR2105872_R1.fastq.gz
│       └── SRR2105872_R2.fastq.gz
├── genome.fasta
```

- `3d_genome_n_crassa` is the name of the working directory as specified in the config file `config.yml`.
- Paired-end fastq files are in the directory `3d_genome_n_crassa/fastq_files/FASTQ_ID` with `FASTQ_ID` the identifier of the paired fastq files. Ftasq identifiers are reported in the config file (`config.yml`).


## Build model

Run 3D model construction:

```bash
snakemake --configfile config.yml --cores 4 --use-singularity --use-conda
```

or with debugging options:

```bash
snakemake --configfile config.yml --cores 4 --use-singularity --use-conda -p --verbose
```

## Example: build model for *Neurospora crassa*

Download and prepare genome sequence:

```bash
bash prepare_n_crassa_genome.sh
```

Define parameters in `config.yml`:

```
workdir: "3d_genome_n_crassa"

organism: "Neurospora crassa"

sra_ids:
- SRR2105869
- SRR2105870
- SRR2105871
- SRR2105872

hicpro_restriction_sites: "dpnii T^TAA"

hicpro_resolutions:
- 10000
- 20000
- 40000
- 50000

pastis_resolutions:
- 50000
- 40000
```

Run the 3D model construction:

```bash
snakemake --configfile config.yml --cores 4 --use-singularity --use-conda
```


## Build DAG graph

For visualisation purpose, you can build the graph of all computational steps involved in the 3D construction of the genome.

```bash
snakemake --configfile config.yml --rulegraph  | dot -Tpdf > rules.pdf
```

With wildcards:

```bash
snakemake --configfile config.yml --dag  | dot -Tpdf > dag.pdf
```


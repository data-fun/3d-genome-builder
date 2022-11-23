<img align="right" width="200px" 
    src="assets/Neurospora_crassa_WT_50kb.gif"
    alt="3D structure of the Neurospora crassa genome at 50 kb resolution">

# 3D genome builder

3d genome builder (3GDB) is a integrated solution to build genome 3D models from HiC raw data and visually integrate omics data on them.
3DGB bundles HiC-Pro and Pastis-nb with additional PDB output file formatting steps into a unified Snakemake workflow with limited input (see *Preparing Required Files*) and unified results in an html output file.

## Download this repository

```bash
git clone https://github.com/data-fun/3d-genome-builder.git
cd 3d-genome-builder
```

## Install dependencies

### Singularity

Download the latest version [here](https://github.com/apptainer/singularity/releases)

Install Singularity:

```bash
sudo apt install -y ./singularity-container_3.8.7_amd64.deb
```

Verify version:

```
$ singularity --version
singularity version 3.8.7
```

### Conda environment

Install [conda](https://docs.conda.io/en/latest/miniconda.html).

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

### Download  HiC-Pro Singularity image


```bash
wget --ciphers=DEFAULT:@SECLEVEL=1 https://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_3.1.0_ubuntu.img -P images
```

Verify HiC-Pro version with:

```bash
$ singularity exec images/hicpro_3.1.0_ubuntu.img HiC-Pro --version
[...]
HiC-Pro version 3.1.0
```

and bowtie2 version:

```bash
$ singularity exec images/hicpro_3.1.0_ubuntu.img bowtie2 --version  2>/dev/null | head -n 1
/usr/local/conda/envs/hicpro/bin/bowtie2-align-s version 2.4.4
```


## Prepare required files

### Create the config file

Create and edit a configuration file in [yaml](https://en.wikipedia.org/wiki/YAML) format. See for instance the template `config_template.yml`

### Add the reference genome

The reference genome fasta file must be located in `WORKING_DIR/genome.fasta` where `WORKING_DIR` is the name of the working directory as specified in your config file.

### Add FASTQ files (optional)

If you already have fastq files stored locally or some fastq files are not available on GEO or SRA, you can use these files providing they are in the proper directory structure:

```
WORKING_DIR/
├── fastq_files
│   ├── ID1
│   │   ├── ID1_R1.fastq.gz
│   │   └── ID1_R2.fastq.gz
│   ├── ID2
│   │   ├── ID2_R1.fastq.gz
│   │   └── ID2_R2.fastq.gz
│   ├── ID3
│   │   ├── ID3_R1.fastq.gz
│   │   └── ID3_R2.fastq.gz
│   └── ID4
│       ├── ID4_R1.fastq.gz
│       └── ID4_R2.fastq.gz
├── genome.fasta
```

- `WORKING_DIR` is the name of the working directory as specified in your config file.
- Paired-end fastq files are in the directory `WORKING_DIR/fastq_files/IDx` with `IDx` the identifier of the paired fastq files. Fastq identifiers are reported in the config file. Please note fastq files have to follow the pattern `<sample ID>_R<1 or 2>.fastq.gz`.

> **Note**
>
> Please strictly follow this file organization as it is later expected by the 3D model construction workflow.

## Build model

Run 3D model construction:

```bash
snakemake --profile smk_profile -j 4 --configfile YOUR-CONFIG.yml
```

> **Note**
> - Adapt `YOUR-CONFIG.yml` to the exact name of the config file you created.
> - Option `-j 4` tells Snakemake to use up to 4 cores. If you are more cores available, you can increase this value (*e.g.* `-j 16`).

Or with debugging options:

```bash
snakemake --profile smk_profile -j 4 --configfile YOUR-CONFIG.yml --verbose
```

Depending on the number and size of fastq files, the 3D construction will take a couple of hours to run.

For troubleshooting, have a look to log files in `WORKING_DIR/logs`, where `WORKING_DIR` is the name of the working directory as specified in your config file.

## Map quantitative values to the 3D model

To add quantitative values to the model run :

```
python ./scripts/map_parameter.py --pdb "path/to/structure.pdb" --BedGraph path/to/annotation.bedgraph --output path/to/output.pdb
```

The quantitative values need to be in a bedgraph file formatted with 4 columns (chromosome/start/stop/value):

```
chr1	0	50000	116.959
chr1	50000	100000	48.4495
chr1	100000	150000	22.8726
chr1	150000	200000	84.3106
chr1	200000	250000	113.109
```

## Interpolate genes into the 3D model

Build a 3D model with one bead for each genes in a given bedgraph :

```
python ./scripts/genes_interpol.py --pdb path/to/structure.pdb --fasta path/to/genome.fasta --resolution HiC-resolution --annotation path/to/genes.bedgraph --output path/to/output.pdb
```

Please adapt `path/to/structure.pdb`, `path/to/genome.fasta`, `HiC-resolution`, `path/to/genes.bedgraph` and `path/to/output.pdb` to your own settings.

The gene list need to be in a bedgraph file formatted with 4 tab-separated columns (chromosome/start/stop/gene_ID):

```
1	1988	1990	NCU10129
1	3386	3388	NCU09901
1	9604	9606	NCU09903
1	15930	15932	NCU11134
1	17872	17874	NCU09904
```

## Examples

- [Build wild type model for *Neurospora crassa*](examples/n_crassa.md)

## Build DAG graph

For visualisation purpose, you can build the graph of all computational steps involved in the 3D construction of the genome.

```bash
snakemake --profile snakemake_profile --rulegraph  | dot -Tpdf > rules.pdf
```

With wildcards:

```bash
snakemake --profile snakemake_profile --dag  | dot -Tpdf > dag.pdf
```


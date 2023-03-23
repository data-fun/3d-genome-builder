<img align="right" width="200px" 
    src="assets/Neurospora_crassa_WT_50kb.gif"
    alt="3D structure of the Neurospora crassa genome at 50 kb resolution">

# 3D genome builder (3DGB)

3D genome builder (3DGB) is a workflow to build 3D models of genomes from HiC raw data and to integrate omics data on the produced models for further visual exploration.
3DGB bundles [HiC-Pro](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x), [PASTIS](https://academic.oup.com/bioinformatics/article/30/12/i26/385087) and custom Python scripts into a unified Snakemake workflow with limited inputs (see *Preparing Required Files*). 3DGB produces annotated 3D models of genome in PDB and G3D formats.

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
└── genome.fasta
```

- `WORKING_DIR` is the name of the working directory as specified in your config file.
- Paired-end fastq files are in the directory `WORKING_DIR/fastq_files/IDx` with `IDx` the identifier of the paired fastq files. Fastq identifiers are reported in the config file. Please note fastq files have to follow the pattern `<sample ID>_R<1 or 2>.fastq.gz`.

> **Note**
>
> Please strictly follow this file organization as it is required by the 3DGB workflow.

## Build model

Run 3DGB:

```bash
snakemake --profile smk_profile -j 4 --configfile YOUR-CONFIG.yml
```

> **Note**
> - Adapt `YOUR-CONFIG.yml` to the exact name of the config file you created.
> - Option `-j 4` tells Snakemake to use up to 4 cores. If you are more cores available, you can increase this value (*e.g.* `-j 16`).

Or with debugging options:

```bash
snakemake --profile smk_profile_debug -j 4 --configfile YOUR-CONFIG.yml --verbose
```

Depending on the number and size of fastq files, the 3D construction will take a couple of hours to run.

For troubleshooting, have a look to log files in `WORKING_DIR/logs`, where `WORKING_DIR` is the name of the working directory as specified in your config file.

## Map quantitative values on the 3D model

To map quantitative values on the model run:

```bash
python ./scripts/map_parameter.py --pdb path/to/structure.pdb --bedgraph path/to/annotation.bedgraph --output path/to/output.pdb
```

Quantitative values should be formatted in a 4-column bedgraph file (chromosome/start/stop/value):

```
chr1	0	50000	116.959
chr1	50000	100000	48.4495
chr1	100000	150000	22.8726
chr1	150000	200000	84.3106
chr1	200000	250000	113.109
```

Each bead of the model will be assigned a quantitative value. The resolution in the bedgraph file should match the resolution used to build the model.


## Get results

Upon completion, the `WORKING_DIR` should look like this:

```
WORKING_DIR/
├── contact_maps
├── dense_matrix
├── fastq_files
├── HiC-Pro
├── logs
├── pastis
├── sequence
└── structure
```

The following paths contain the most interesting results:

- `WORKING_DIR/contact_maps/*.png` : contact maps.
- `WORKING_DIR/HiC-Pro/output/hic_results/pic/*/*.pdf` : graphical summaries of read alignments produced by Hi-C Pro.
- `WORKING_DIR/pastis/structure_RESOLUTION.pdb` : raw 3D models (in PDB format) produced by Pastis.
- `WORKING_DIR/structure/RESOLUTION/structure_cleaned.*` : final (annotated) 3D models in PDB and G3D formats.

> **Note**
> - `WORKING_DIR` is the name of the working directory as specified in your config file.
> - `RESOLUTION` is the resolution of the Hi-C data specified in the config file.

## Examples

- [Wild type model for *Neurospora crassa*](examples/n_crassa.md)
- [Models built for the 3DGB paper](examples/paper/paper.md)

## Visualize 3D model structures

To visualize 3D model structures (.pdb and .g3d files), follow this quick [tutorial](visualization/visualization.md).


## Build DAG graph

For visualization purpose, you can build the graph of all computational steps involved in the 3D construction of the genome.

```bash
snakemake --profile smk_profile --configfile YOUR-CONFIG.yml --rulegraph  | dot -Tpdf > rules.pdf
```

where `YOUR-CONFIG.yml` should be replaced by the name of the config file you created.

With wildcards:

```bash
snakemake --profile smk_profile --configfile YOUR-CONFIG.yml --dag  | dot -Tpdf > dag.pdf
```


# 3D genome builder

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

## Download  HiC-Pro Singularity image


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


## Add parameters

Create and edit the configuration file `config.yml` in yaml format. See for instance `config_template.yml`

## Add reference genome

The reference genome fasta file must be located in `WORKING_DIR/genome.fasta` where `WORKING_DIR` is the name of the working directory as specified in the config file `config.yml`.

## Add FASTQ files (optional)

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

- `WORKING_DIR` is the name of the working directory as specified in the config file `config.yml`.
- Paired-end fastq files are in the directory `WORKING_DIR/fastq_files/IDx` with `IDx` the identifier of the paired fastq files. Fastq identifiers are reported in the config file (`config.yml`). Please note fastq files have to follow the pattern `<sample ID>_R<1 or 2>.fastq.gz`.

> **Note**
>
> Please strictly follow this file organization as it is later expected by the 3D model construction workflow.

## Build model

Run 3D model construction:

```bash
snakemake --profile snakemake_profile -j 4
```

> **Note**
>
> Option `-j 4` tells Snakemake to use up to 4 cores. If you are more cores available, you can increase this value (*e.g.* `-j 16`).

Or with debugging options:

```bash
snakemake --profile snakemake_profile -j 4 --verbose
```

Depending on the number and size of fastq files, the 3D construction will take a couple of hours to run.

## Example: build model for *Neurospora crassa*

1. Download and prepare the reference genome sequence:

    ```bash
    bash prepare_n_crassa_genome.sh
    ```

2. Copy parameters from [config_n_crassa.yml](config_n_crassa.yml):

    ```bash
    cp config_n_crassa.yml config.yml
    ```

3. Run the 3D model construction:

    ```bash
    snakemake --profile snakemake_profile -j 4
    ```

## Build DAG graph

For visualisation purpose, you can build the graph of all computational steps involved in the 3D construction of the genome.

```bash
snakemake --profile snakemake_profile --rulegraph  | dot -Tpdf > rules.pdf
```

With wildcards:

```bash
snakemake --profile snakemake_profile --dag  | dot -Tpdf > dag.pdf
```


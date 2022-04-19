# 3D genome builder

## Download this repository

```bash
$ git clone https://github.com/data-fun/3d-genome-builder.git
$ cd 3d-genome-builder
```

## Install dependencies

### Singularity

Download last version [here](https://github.com/apptainer/singularity/releases)

Install package:

```bash
$ sudo apt install -y ./singularity-container_3.8.7_amd64.deb
$ singularity --version
singularity version 3.8.7
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

## Build HiC-Pro Singularity image

```bash
git clone https://github.com/nservant/HiC-Pro
cd HiC-Pro
sudo singularity build ../images/hicpro.sif Singularity
cd ..
rm -rf HiC-Pro
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

Create and edit the configuration file (`config.yml`).


## Build model

Download and prepare genome sequence:

```bash
bash prepare_n_crassa_genome.sh
```

Run 3D model construction:

```bash
snakemake --cores 4 --configfile config.yml
```

or for debugging purpose:

```bash
snakemake --cores 4 --configfile config.yml -p --verbose --debug-dag
```

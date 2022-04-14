# 3D genome builder


## Install dependencies

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
sudo singularity build ../images/hicpro.img Singularity
cd ..
rm -rf HiC-Pro
```

## Add parameters

Create and edit the configuration file (`config.yml`).

## Build mode

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

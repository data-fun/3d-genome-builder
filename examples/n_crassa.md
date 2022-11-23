# Build a 3D genome model for *Neurospora crassa*


## Download and prepare the reference genome sequence


```bash
bash examples/n_crassa_prepare_genome.sh
```

## Build 3D genome model

```bash
snakemake --profile smk_profile -j 4 --configfile examples/n_crassa_WT.yml
```
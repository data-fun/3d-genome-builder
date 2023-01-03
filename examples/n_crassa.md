# Build a 3D genome model for *Neurospora crassa*

The HiC data are from [Galazka et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27260477/) (SRA ID : [SRX1099807](https://www.ncbi.nlm.nih.gov/sra?term=SRX1099807)) and the ChipSeq data from [Jamieson et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26537359/) (SRA ID : [SRR2026390](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2026390)) and [Basenko et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26578794/) (SRA ID : [SRR2036168](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2036168)).


## Download and prepare the reference genome sequence

```bash
bash examples/n_crassa_prepare_genome.sh
```

## Build 3D genome model

```bash
snakemake --profile smk_profile -j 4 --configfile examples/n_crassa_WT.yml
```

See output file `3DGB_n_crassa_WT/structure/50000/structure_cleaned.pdb` that should be equivalent to `./examples/structure_cleaned.pdb`.

## Map ChipSeq values to the 3D model

```bash
python scripts/map_parameter.py --pdb "./3DGB_n_crassa_WT/structure/50000/structure_cleaned.pdb" --BedGraph "./examples/n_crassa.bedgraph" --output "./3DGB_n_crassa_WT/structure/50000/structure_with_parameter.pdb"
```

See output file `./3DGB_n_crassa_WT/structure/50000/structure_with_parameter.pdb` that should be equivalent to `./examples/structure_with_parameter.pdb`.

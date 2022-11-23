# Build a 3D genome model for *Neurospora crassa*

HiC data are from [Galaska et al, 2016](https://pubmed.ncbi.nlm.nih.gov/27260477/) (SRA experiment [SRX1099807](https://www.ncbi.nlm.nih.gov/sra?term=SRX1099807)).

ChIP-Seq data (for mapping) are from [Jamieson et al, 2016](https://pubmed.ncbi.nlm.nih.gov/26537359/) (SRA run [SRR2026390](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2026390)) and [Basenko et al, 2015](https://pubmed.ncbi.nlm.nih.gov/26578794/) (SRA run [SRR2036168](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2036168)).

## Download and prepare the reference genome sequence

```bash
bash examples/n_crassa_prepare_genome.sh
```

## Build 3D genome model

```bash
snakemake --profile smk_profile -j 4 --configfile examples/n_crassa_WT.yml
```

## Map ChIP-Seq values to the 3D model

```bash
python scripts/map_parameter.py --pdb "./3DGB_n_crassa_WT/structure/50000/structure_cleaned.pdb" --BedGraph "./examples/n_crassa.bedgraph" --output "./3DGB_n_crassa_WT/structure/50000/structure_with_parameter.pdb"
```

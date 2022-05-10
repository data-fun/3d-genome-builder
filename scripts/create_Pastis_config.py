from pathlib import Path

import jinja2
import pandas as pd


# Create genome size with chromosome size, without chromosome number
chromosome_sizes_df = pd.read_csv(snakemake.input.chromosome_sizes, 
                                  sep="\t", header=None, names=["chromosome", "size"])
chromosome_number = chromosome_sizes_df.shape[0]
chromosome_sizes_df["size"].to_csv(snakemake.output.chromosome_sizes, header=False, index=False)

# Create Pastis config file
with open(snakemake.input.template, "r") as template_file, \
     open(snakemake.output.config, "w") as out_file:
    template = jinja2.Template(template_file.read())
    out_file.write(template.render(
        resolution=snakemake.wildcards.resolution,
        chromosomes=",".join([str(chrom_id+1) for chrom_id in range(chromosome_number)]),
        sizes=Path(snakemake.output.chromosome_sizes).name,
        npy=Path(snakemake.input.npy).resolve()
        ))

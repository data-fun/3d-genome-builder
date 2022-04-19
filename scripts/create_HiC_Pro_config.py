from pathlib import Path

import jinja2

with open(snakemake.input[0], "r") as template_file, \
     open(snakemake.output[0], "w") as out_file:
    template = jinja2.Template(template_file.read())
    out_file.write(template.render(
        genome_index_path=Path(snakemake.params.genome_index_path).resolve(),
        chromosome_sizes=Path(snakemake.input.chromosome_sizes).resolve(),
        genome_fragment=Path(snakemake.input.genome_fragment).resolve()
        ))

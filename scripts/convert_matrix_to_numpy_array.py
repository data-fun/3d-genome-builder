import numpy as np

raw_matrix = np.loadtxt(snakemake.input[0])
np.save(snakemake.output[0], raw_matrix)

import matplotlib.pyplot as plt
import numpy as np


contact_map = np.loadtxt(snakemake.input[0])

plt.imshow(np.log2(contact_map+1), cmap='hot')
plt.colorbar()

plt.savefig(snakemake.output[0])
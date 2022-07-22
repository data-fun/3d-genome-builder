# Post-deployment script
# https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

# Clone pastis repo
git clone --depth 1 https://github.com/hiclib/pastis.git
# Copy required submodules from the last commit
cp -R pastis/pastis/dispersion "${CONDA_PREFIX}/lib/python3.9/site-packages/pastis/"
cp -R pastis/pastis/optimization "${CONDA_PREFIX}/lib/python3.9/site-packages/pastis/"
# Remove pastis repo
rm -rf pastis

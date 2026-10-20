# partis with R and TreeSim, for 'partis simulate --simulate-from-scratch'.
#
# partis draws each clonal family's tree with the R package TreeSim
# (treegenerator.run_treesim through utils.run_r), and the partis image has
# no R (CI run 6: 'command R not found in path'). This image adds R and
# TreeSim in their own micromamba environment, so the partis environment is
# not changed, and puts that environment last on PATH, where partis finds R.
#
# The base is the partis image by digest from images.lock (build argument
# PARTIS_BASE, set by 'make pin'). R and TreeSim come from conda-forge only,
# fixed by version and build string. The full solved package list (name,
# version, build, channel) is kept in the image at
# /opt/conda/envs/r/packages.txt, and
# 'make pin' records the built image's ID in images.lock.
ARG PARTIS_BASE
FROM ${PARTIS_BASE}
ARG R_BASE_SPEC=r-base=4.5.3=h502d0c9_4
ARG R_TREESIM_SPEC=r-treesim=2.4=r45hc72bb7e_5
RUN command -v micromamba \
 && micromamba create -y -n r --override-channels -c conda-forge "${R_BASE_SPEC}" "${R_TREESIM_SPEC}" \
 && micromamba list -n r > /opt/conda/envs/r/packages.txt \
 && micromamba clean --all --yes
ENV PATH=$PATH:/opt/conda/envs/r/bin
# The build stops here if R does not start or TreeSim does not load.
RUN R --version \
 && Rscript -e 'library(TreeSim); cat("TreeSim", format(packageVersion("TreeSim")), "\n")'

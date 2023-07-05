FROM continuumio/miniconda3:latest

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding anaconda to the list of channels
RUN conda config --add channels anaconda

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Adding conda-forge to the list of channels
RUN conda config --add channels conda-forge

# Installing mamba
RUN conda install -y mamba

# Installing packages
RUN mamba install -y \
    h5py \
    python-igraph \
    matplotlib \
    numpy \
    pandas \
    scanorama \
    scipy \
    scikit-learn \
    seaborn \
    statsmodels \
    umap-learn && \
    conda clean -afty

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]

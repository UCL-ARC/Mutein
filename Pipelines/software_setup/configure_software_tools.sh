#!/usr/bin/env bash

#
# setup conda environments for the pipeline
#

set -eu
source ~/.mutein_settings

#make conda commands available
module load ${MUT_CONDA_MODULE}
echo ". ${MUT_CONDA_SETUP}" >> ~/.bashrc

#add channels contains bioinformatics tools available
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#required for enaDataGet command
conda create --name enabrowsertools
conda activate enabrowsertools
conda install enabrowsertools
conda deactivate

#required for trim_galore
conda create --name trim-galore
conda activate trim-galore
conda install trim-galore
conda deactivate

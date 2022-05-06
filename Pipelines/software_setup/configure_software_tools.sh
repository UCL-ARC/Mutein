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
conda create --yes --name enabrowsertools
conda activate enabrowsertools
conda install --yes "enabrowsertools>=1.5.4"
conda deactivate

#required for trim_galore
conda create --yes --name trim-galore
conda activate trim-galore
conda install --yes "trim-galore>=0.6.7"
conda deactivate

#required for bwa
conda create --yes --name bwa
conda activate bwa
conda install --yes "bwa>=0.7.15"
conda install --yes "samtools>=1.15.1"
conda deactivate

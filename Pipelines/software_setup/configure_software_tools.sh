#!/usr/bin/env bash

#
# setup conda environments for the pipeline
# before running this script ensure you have installed and configured miniconda
# and have copied example_mutein_settings to ~/.mutein_settings
# and edited it appropriately

# typically to install miniconda you will download and run the installer:
#   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#   chmod u+x ./Miniconda3-latest-Linux-x86_64.sh
#   ./Miniconda3-latest-Linux-x86_64.sh
# then log out and in again after miniconda has inserted itself into your bashrc setup
# then copy and edit example_mutein_settings to ~/.mutein_settings before proceeding

set -eu
source ~/.mutein_settings

#add channels contains bioinformatics tools available
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#required for enaDataGet command
conda create --yes --name ${MUT_PREFIX}enabrowsertools
conda activate ${MUT_PREFIX}enabrowsertools
conda install --yes "enabrowsertools>=1.5.4"
conda deactivate

#required for trim_galore
conda create --yes --name ${MUT_PREFIX}trim-galore
conda activate ${MUT_PREFIX}trim-galore
conda install --yes "trim-galore>=0.6.7"
conda deactivate

#required for bwa
conda create --yes --name ${MUT_PREFIX}bwa
conda activate ${MUT_PREFIX}bwa
conda install --yes "bwa>=0.7.15"
conda install --yes "samtools>=1.15.1"
conda deactivate

#required for gatk
conda create --yes --name ${MUT_PREFIX}gatk4
conda activate ${MUT_PREFIX}gatk4
conda install --yes gatk4
conda deactivate

#required for pyega3
conda create --yes --name ${MUT_PREFIX}pyega3
conda activate ${MUT_PREFIX}pyega3
conda install --yes pyega3
conda deactivate

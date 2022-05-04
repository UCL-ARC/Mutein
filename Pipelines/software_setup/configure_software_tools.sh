#!/usr/bin/env bash

#
# setup conda environments for the pipeline
#

#source mutein settings
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source ${SCRIPT_DIR}/../software_setup/mutein_settings.sh

#load miniconda as a module or setup miniconda manually if no module exists
module load python/miniconda3/${MUTEIN_CONDA_VER}

#one off conda config to set it up fully for your user account
echo ". /shared/ucl/apps/miniconda/${MUTEIN_CONDA_VER}/etc/profile.d/conda.sh" >> ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#required for enaDataGet
conda create --name enabrowsertools
conda activate enabrowsertools
conda install enabrowsertools
conda deactivate

#required for trim_galore
conda create --name trim-galore
conda activate trim-galore
conda install trim-galore
conda deactivate

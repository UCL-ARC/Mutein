#!/bin/bash

#export PATH=${PATH}:/opt/singularity/bin
export DEFFILE=/home/ccaervi/Mutein/Pipelines/singularity/mutein.def

sudo /opt/singularity/bin/singularity build --allow-setuid mutein.sif ${DEFFILE}

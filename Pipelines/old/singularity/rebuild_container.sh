#!/bin/bash

#
# rebuild the mutein.sif container image
# needs to be run on a linux machine where you have root
# eg an ubuntu linux vm in virtualbox on your local machine
# you will need to edit the paths to the sif and def file and singularity binary
#

#export PATH=${PATH}:/opt/singularity/bin
export DEFFILE=/home/ccaervi/Mutein/Pipelines/singularity/mutein.def
export SIFFILE=/home/ccaervi/549_mutein/containers/mutein.sif
export BINFILE=/opt/singularity/bin/singularity

sudo ${BINFILE} build ${SIFFILE} ${DEFFILE}

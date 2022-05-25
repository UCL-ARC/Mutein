#!/bin/bash -l
module load python3/recommended

mode=$1
pattern=$2
WorkDir=$3			#/home/ucbtlcr/Scratch/workspace/
DataDir=$4			#/home/ucbtlcr/MuteinData/
InstallDir=$5		#/home/ucbtlcr/Mutein/
PipelineDir=$6		#/home/ucbtlcr/Mutein/Pipelines/geneanalysis/

script=${InstallDir}Pipelines/libs/remote.py
python ${script} $mode $pattern $WorkDir $DataDir $InstallDir $PipelineDir


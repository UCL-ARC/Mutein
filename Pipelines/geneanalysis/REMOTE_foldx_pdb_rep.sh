#!/bin/bash -l
module load python3/recommended

run=qsub

dataset=$1
gene=$2
pdb=$3
WorkDir=$4			#/home/ucbtlcr/Scratch/workspace/
DataDir=$5			#/home/ucbtlcr/MuteinData/
InstallDir=$6		#/home/ucbtlcr/Mutein/
PipelineDir=$7		#/home/ucbtlcr/Mutein/Pipelines/geneanalysis/

config=${PipelineDir}config/batch_gene_1_rep.yml
script=${InstallDir}Pipelines/libs/pipeline_qsubber.py

cd $DataDir
echo "Submitting ${script} $install_dir $PWD ${config} $run $dataset"
#install_dir, working_dir, yaml_file, py_or_sh, dataset, gene, pdb
python ${script} $InstallDir $DataDir ${config} $run $dataset "" ""


 
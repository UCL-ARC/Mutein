#!/bin/bash -l
module load python3/recommended

jobid=$1
install_dir=$2
script=${install_dir}Pipelines/foldx/libs/pipeline_view.py

echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"

echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} "JOBID=$jobid"


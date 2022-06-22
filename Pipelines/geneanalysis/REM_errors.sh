#!/bin/bash -l
module load python3/recommended

mode=$1
install_dir=$2
script=${install_dir}Pipelines/libs/pipeline_delete.py

echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"

echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} "MODE=$mode"




# echo ""
# echo ""
# echo "Current directory contents"
# target="/home/ucbtlcr/Scratch/workspace/"
# let count=0
# for f in "$target"/*
# do
#     let count=count+1
# done
# echo "Count: $count"

# outputString=$(python ${script} "MODE=$mode")
# echo $outputString

# target="/home/ucbtlcr/Scratch/workspace/"
# let count=0
# for f in "$target"/*
# do
 #    echo $(basename $f)
  #   let count=count+1
# done
# echo ""
# echo "Count: $count"


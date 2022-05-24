#!/bin/bash

### can be clean or all

module ()
{
    eval `/shared/ucl/apps/modules/3.2.6/Modules/$MODULE_VERSION/bin/modulecmd bash $*`
}

module load python3/recommended

target="/home/ucbtlcr/Scratch/workspace/"
let count=0
for f in "$target"/*
do
    let count=count+1
done
echo "Count: $count"

outputString=$(python ${script} "MODE=$mode")
echo $outputString

target="/home/ucbtlcr/Scratch/workspace/"
let count=0
for f in "$target"/*
do
    echo $(basename $f)
    let count=count+1
done
echo ""
echo "Count: $count"


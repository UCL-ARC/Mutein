### can be clean or all
mode=$1
install_dir=$2
script=${install_dir}Pipelines/libs/pipeline_delete.py

# assume we are in the home directory and the data is in MuteinData

module load python3/recommended

cd MuteinData

#echo "EXE PATH=$install_dir"
#echo "CURRENT=$PWD"
#echo "SCRIPT=$script"

#echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

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


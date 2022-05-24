### can be clean or all
mode=$1
install_dir=$2
script=${install_dir}Pipelines/libs/pipeline_delete.py

# assume we are in the home directory and the data is in MuteinData
cd MuteinData

echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"
echo "SCRIPT=$script"

echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} "MODE=$mode"
echo "MUTEIN SCRIPT ENDED"


target="/home/ucbtlcr/Scratch/workspace/"
let count=0
for f in "$target"/*
do
    echo $(basename $f)
    let count=count+1
done
echo ""
echo "Count: $count"


### SCRIPT AUTOMATICALLY GENERATED BY MUTEIN PIPELINE (UCL-ARC 2022) ###
### can be clean or all
mode=$1
install_dir=$2
script=${install_dir}Pipelines/libs/pipeline_delete.py

echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"

echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} "MODE=$mode"
echo "MUTEIN SCRIPT ENDED"

### SCRIPT FOT MUTEIN PIPELINE (UCL-ARC 2022) ###
run=qsub
dataset=$1
gene=$2
install_dir=$3
script=${install_dir}Pipelines/libs/pipeline_qsubber.py
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_4_tasks.yml
echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"

echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} $install_dir $PWD ${config} $run $dataset $gene ""




### SCRIPT AUTOMATICALLY GENERATED BY MUTEIN PIPELINE (UCL-ARC 2022) ###
run=qsub
install_dir=$1
script=${install_dir}Pipelines/libs/pipeline_qsubber.py
config=${install_dir}Pipelines/geneanalysis/config/batch_pdb.yml
echo "EXE PATH=$install_dir"
echo "CURRENT=$PWD"
module load python3/recommended
echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python ${script} $install_dir $PWD ${config} $run shearwater FAT1 smhom_6vg4.1.a_2923_3537
python ${script} $install_dir $PWD ${config} $run shearwater FAT1 smhom_6vg1.1.a_2923_3535
python ${script} $install_dir $PWD ${config} $run shearwater FAT1 smhom_6vg4.1.a_2862_3432
python ${script} $install_dir $PWD ${config} $run shearwater FAT1 smhom_6vg4.1.a_767_1346

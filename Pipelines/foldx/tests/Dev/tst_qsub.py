"""
RSA: 19/5/22
This CI test script runs from the pdb level
It enables debugging of the scripts as if run from a batch

"""
import sys, os
# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
config_path = "/".join(dirs) + "/config"
python_path = "/".join(dirs) + "/python"
lib_path = "/".join(dirs) + "/libs"
sys.path.append(python_path)
sys.path.append(lib_path)
sys.path.append(config_path)
import Paths
######################################################################

script=lib_path + "/pipeline_qsubber.py"
data_dir="/home/rachel/UCL/github/MuteinData/"
install_dir="/home/rachel/UCL/github/Mutein/"

config=config_path + "/batch_gene_4_tasks.yml"
#config=config_path + "/batch_gene_4_untasks.yml"
#config=config_path + "/batch_gene_3_prep.yml"
#config=config_path + "/batch_dataset_2_stitch.yml"

dataset = "one"
gene = "brca2"
pdb=""

import pipeline_qsubber as pq
#install_dir, working_dir, yaml_file, py_or_sh, dataset, gene, pdb
pq.pipeline_qsubber(["",install_dir,data_dir,config,"py",dataset,gene,pdb])




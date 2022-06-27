#include file of variables requiring python code rather than just plain YAML
import os
import glob

## ======== conda environments ========    
#generate list of required conda environments from the list of "env" files
conda_env_path = os.path.join(os.environ['MUT_DIR'],"Pipelines/conda_envs")
conda_envs = glob_wildcards(os.path.join(conda_env_path,"{envname}.env"))
all_conda_envs = [f"config/conda_envs/{envname}.touch" for envname in conda_envs.envname]

## ======== reference genome ========
#define the main reference genome fasta path and all its index files
bgz_fasta = os.path.join(config["ref"]["dir"],config["ref"]["fasta"])
exts = ['gzi','fai','dict','bwt','pac','ann','amb','sa']
all_ref_files = [ bgz_fasta ] + list(expand('{base}.{ext}',base=bgz_fasta,ext=exts))

#define just the bwa related index files
bwt_indexes = list(expand('{base}.{ext}',base=bgz_fasta,ext=['bwt','pac','ann','amb','sa'])) 

## ======== datasets ========
#define names of subfolders of the main dataset folder, one per dataset
all_datasets = \
[
    "keogh2018",
    #"yokoyama2019",
]

#path to a "touch" file used to flag that data download is completed
#used because we don't know the names of all the datafiles in advance
#therefore we cannot list them for snakemake
all_datasets_downloaded = [ f"datasets/{x}/download_success" for x in all_datasets ]
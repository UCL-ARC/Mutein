# Run Pipeline on Myriad


#### 1) clone the git repo Mutein, or navigate to it and do git pull, potentially to a branch
```
cd Mutein
git pull
git checkout rachel/foldxv1
```
#### 2) Navigate to the Pipelines directory and run the script, first installing python mods
##### There are 3 choices:
####### 1) The script itself
####### 2) The config
####### 3) "py" to run as a serial python script, "qsub" to submit qsub batches for parallel
```
module load python3/recommended
python3 overall_rsa.py batch_tst01.yml py
```
#### 3) GeneProt: This starts with I think a vcf file and ends with the protein structures
#### 4) FoldX: This pipeline takes a pdb file and a list of mutations and creates xxxx
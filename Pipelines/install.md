# Run Pipeline on Myriad


#### 1) clone the git repo Mutein, or navigate to it and do git pull, potentially to a branch
```
cd Mutein
git pull
git checkout rachel/foldxv1
```
#### 2) Navigate to the Pipelines directory and run the script, first installing python mods
```
module load python3/recommended
```
##### There are a few choices:
```
sh pipeline_covid.sh
sh pipeline_notch.sh
sh pipeline_notch_NOTCH1.sh
```
#### 3) GeneProt: This starts with I think a vcf file and ends with the protein structures
#### 4) FoldX: This pipeline takes a pdb file and a list of mutations and creates xxxx
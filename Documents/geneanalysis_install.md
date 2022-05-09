##### Owner:Rachel
--------------------------------------------------------------------------
## HPC FoldX Pipeline for myriad
--------------------------------------------------------------------------
- When everything is installed and running, the batch scripts are in the main scripts directory as symlinks
- Navigate to your chosen working directory, from wich the inputs/outputs directory are assumed to be present
- Then run the approriate batch symlink
```
>> module load python3/recommended
>> cd ~/MuteinData
>> /home/ucbtlcr/Mutein/Pipelines/geneanalysis/ppl_notch_NOTCH1.sh /home/ucbtlcr/Mutein/

(when the environment variables are created it will be simplified to)
>> ppl_notch_NOTCH1.sh
```
--------------------------------------------------------------------------
### 1.The data is structured in folders
```
scripts/*
inputs/pdb/*
thruputs/pdb/*
results/setname/*
outputs/setname/*
(The thruputs, results and outputs are not synced with github)
```
### 2.Pipeline inputs
##### A pdb file
##### A list of mutations in a file of the form
```
Mutation,Variant
NA501Y,Alpha
AA570D,Alpha
TA716I,Alpha
SA982A,Alpha
DA1118H,Alpha
DA614G,Alpha
```
--------------------------------------------------------------------------
### 3.To run the pipeline on myriad:
##### clone the repo
```
# Clone the main branch or a branch
>> cd Mutein
# Get the code
>> git checkout .
>> git pull

# The data directory needs to be set up with inputs
(this could be a different githun, for now I have a git_sync directory)

```
-----------------------------------------------------------------------
```
The scripts dependency is:
                                    [Per Gene]:01_Gene_To_Protein
                                              |
                                    [Per PDB]:2:FOLDX REPAIR
                                        |                |
                            [Per PDB]:03a:SPLIT         03b:VARIANT SPLIT
                                        |                |
            [ARRAY JOBS]    [Per PDB]:04a:POS_SCAN       04b:SINGLE_SCAN
                                        |                |
                            [Per PDB]:05a:AGGREGATE     05b:VARIANT AGGREGATE
                                                         |
                            [PDBs->GENE]:SCRIPT 06: Structures->Gene under selection 
 ```
-----------------------------------------------------------------------






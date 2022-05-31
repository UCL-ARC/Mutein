##### Owner:Rachel
-----------------------------------------------------------------------
## HPC FoldX Pipeline for myriad
-----------------------------------------------------------------------
```
The scripts structure is:

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
### Overview
This pipeline compares the background protein folding stability of a given pdb structure against known variants.

- The background folding stability is calculated by mutating every single possible residue in a protein structure
- The variants are calculating by permuting the possibilities of the known variants
-----------------------------------------------------------------------
### How to run
- Navigate to your chosen data directory
- Load python
```
module load python3/recommended
```
- Run the scripts - in this example the notch1 test data
- - they have to be run one by one, wait until they finish. The 1st and 3rd take moments, the 2nd takes.... a day...?
- there are 3 scripts to run in the pipeline
- - notch turns the entire vcs dataset into genes
- - - (notch could automatically kick off all the sub genes in the future)
- - notch_NOTCH1 runs the pdb files for the chosen gene
- - notch_NOTCH1_stitch aggregates the separate structures together
- - - (only not a dependent script due to multiple dependencies at the moment)

### TEST DATA SETS 
# Amend this file for testing
```
/home/ucbtlcr/Mutein/Pipelines/geneanalysis/foldx__scratch.sh /home/ucbtlcr/Mutein/
```
### This is a small dataset 
```
python ${script} $install_dir $PWD ${config} $run "shearwater" "NOTCH1" 1pb5
```
### If the gene is ALL 
- Then it loops through all the genes in a dataset and kicks of batches for every one (this could be quite a lot!!)
- To run this in a sensible order, choose the yaml files that are appropriate
- - DATASET first prepare the data
```
config=${install_dir}Pipelines/geneanalysis/config/batch_dataset_prep.yml
```
- - GENE only
```
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_1_pdbs.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_2_rep.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_3_prep.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_4_tasks.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_gene_5_agg.yml
```
- - PDB only
```
config=${install_dir}Pipelines/geneanalysis/config/batch_pdb_1_rep.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_pdb_2_prep.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_pdb_3_tasks.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_pdb_4_agg.yml
config=${install_dir}Pipelines/geneanalysis/config/batch_dataset_prep.yml
 ```
#### Clean or rerun up at the end to remove successful logs leaving only errors
To control errors and log files this in-process script can be run at any time.
- Identify errors and remove finished happy log files
-- /home/ucbtlcr/Mutein/Pipelines/geneanalysis/ppl_cleanup_and_report.sh CLEAN /home/ucbtlcr/Mutein/
- Rerun failed batches
-- /home/ucbtlcr/Mutein/Pipelines/geneanalysis/ppl_cleanup_and_report.sh ALL /home/ucbtlcr/Mutein/
- Rerun all batches (this would be if evrythong just got disconnected or something, don;t do this if they are still running!)
-- /home/ucbtlcr/Mutein/Pipelines/geneanalysis/ppl_cleanup_and_report.sh RERUNERROR /home/ucbtlcr/Mutein/
- Delete all the log files, no checking, only if you have finished looking at them
-- /home/ucbtlcr/Mutein/Pipelines/geneanalysis/ppl_cleanup_and_report.sh RERUNALL /home/ucbtlcr/Mutein/

### Inputs
1. A pdb file
2. A variants file
-----------------------------------------------------------------------
### Pipeline parent script
- The 00 script is a parent script that runs all the other scripts with the correct dependencies
-----------------------------------------------------------------------
### Scripts
##### prepare
- 1] This takes the gene and finds the appropriate pdb structues
- 2] This runs foldx repairs on a pdb file, the repairs ensure that the atom positions are relaxed into a favoured position.
##### background mutations
- 3a] This splits the possible mutations up into a config file ready for parallelizing
- 4a] An array job: This calls foldx positionscan function with 1 of the rows from the config file, 1 set of mutations
- 5a] This aggregates and performs some analysis on all the individual mutations
##### variant mutations
- 3b] This makes a config for the variants
- 4b] An array job: This calls posscan for the variants
- 5b] And finally aggregates and analyses
##### aggregate to the gene
- 6] Aggregates all the pdb results into a gene result directory

-----------------------------------------------------------------------
### Structure
There are pairs of scripts bash/python. The bash submits to qsub, the python runs in python. This enables a test environment to run the python scripts via both CI and user testing from the tests directory.

There are 3 levels of parameters. In order of least preferred to most:
- Batch params (eg number of array jobs, time allowed, script dependence)
- Pdb parameters (eg override of arrays for variants, name of variant, number of repairs)
- Command line inputs - pdb needed, override anything else if you want (eg change the number of array jobs for testing)
-----------------------------------------------------------------------
### Environment
The scripts can be run in 3 modes - python, hpc and inputs (inputs only records the inputs to external libraries). You can override this with user= in the batch.
You can add your dev environment to the top of helper.py
An example of a batch command with parameters overridden:


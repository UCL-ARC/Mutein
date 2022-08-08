#### Owner:Rachel
-----------------------------------------------------------------------
## HPC FoldX Pipeline for myriad - the scripts to run

-----------------------------------------------------------------------
## Overview of how to run
The scripts can be run in python or qsub mode. Tbhat setting is within the bash script as run=qsub (or run=py)

(n.b. the support gui has a button for each of these scripts).

-------------------------------------------------------------
# Full dataset
Using an example of mouse, and assuming each script is kicked off individually from a linux ssh shell - 

- My install directory is /home/ucbtlcr/Scratch/Mutein/
- My script directory is /home/ucbtlcr/Scratch/Mutein/Pipelines/foldx/scripts/
- My data directory is /home/ucbtlcr/Scratch/MuteinData/Foldx/
- My log directory is /home/ucbtlcr/Scratch/workspace/
- My pipeline directory is /home/ucbtlcr/Scratch/Mutein/Pipelines/foldx/ 

the commands use the directory locations, I will simplify to make it clearer, so:

```
{scripts_dir}REMOTE_foldx_pdb_rep.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
means:
```
/home/ucbtlcr/Scratch/Mutein/Pipelines/foldx/scripts/REMOTE_foldx_pdb_rep.sh mouse Abcb11 smhom_6lr0_1_a_1_1315 /home/ucbtlcr/Scratch/workspace/ /home/ucbtlcr/Scratch/MuteinData/Foldx/ /home/ucbtlcr/Scratch/Mutein/ /home/ucbtlcr/Scratch/Mutein/Pipelines/foldx/ 
```

1. Download genes and pdbs
```
{scripts_dir}REMOTE_foldx_dataset_pdbs.sh DATASET x x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
2. Repair all pdbs
```
{scripts_dir}REMOTE_foldx_gene_rep.sh DATASET ALL x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
3. Create splits file for parallelizing
```
{scripts_dir}REMOTE_foldx_gene_split.sh DATASET ALL x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
4. Submit all tasks
```
{scripts_dir}REMOTE_foldx_gene_tasks.sh DATASET ALL x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
5. Submit all aggregation
```
{scripts_dir}REMOTE_foldx_dataset_agg.sh DATASET x x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```

## utility commands
- view the status of the batches
```
{scripts_dir}REMOTE.sh GENES DATASET:X:X {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
- delete non-error finished log files
```
{scripts_dir}REM_errors.sh CLEAN {install_dir}
```
- There is a bug here as it assumes the position of the log directory as:
    - scratch_dir = "/home/" + homeuser + "/Scratch/workspace/"

# Single gene

You may or may not have a dataset, if not leave it blank.
The presence of a dataset pushes the directory structure down a level.

1. Download genes and pdbs
```
{scripts_dir}REMOTE_foldx_gene_pdbs.sh DATASET GENE x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
2. Repair all pdbs
```
{scripts_dir}REMOTE_foldx_gene_rep.sh DATASET GENE x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
3. Create splits file for parallelizing
```
{scripts_dir}REMOTE_foldx_gene_split.sh DATASET GENE x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
4. Submit all tasks
```
{scripts_dir}REMOTE_foldx_gene_tasks.sh DATASET GENE x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
5. Submit all aggregation
```
{scripts_dir}REMOTE_foldx_gene_agg.sh DATASET GENE x {log_dir} {data_dir} {install_dir} {pipeline_dir}
```

## utility commands
- view the status of the batches
```
{scripts_dir}REMOTE.sh PDBS DATASET:GENE:X {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
# Single pdb structure

You may or may not have a dataset and gene, if not leave them blank.
The presence of a dataset/gene pushes the directory structure down a level.

1. Download pdb
```
{scripts_dir}REMOTE_foldx_pdb_pdb.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
2. Repair pdb
```
{scripts_dir}REMOTE_foldx_pdb_rep.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
3. Create splits file for parallelizing
```
{scripts_dir}REMOTE_foldx_pdb_split.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
4. Submit all tasks
```
{scripts_dir}REMOTE_foldx_pdb_tasks.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```
5. Submit all aggregation
```
{scripts_dir}REMOTE_foldx_pdb_agg.sh DATASET GENE PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```

## utility commands
- view the status of the batches
```
{scripts_dir}REMOTE.sh PDB DATASET:GENE:PDB {log_dir} {data_dir} {install_dir} {pipeline_dir}
```

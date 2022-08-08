##### Owner:Rachel
-----------------------------------------------------------------------
## VCF to Structure Pipeline
-----------------------------------------------------------------------
##### This pipeline is not yet complete.

Note that as yet there has not been a need to make any of these into array jobs (parallelised).
Although they can take some time (a couple of hours) that is not longer than an optimum array job.
```
The scripts structure is:
(* are not yet done)

            *SCRIPT 01: VCF->Genes under selection (NOT DONE)
                                |               
            SCRIPT 02: GENES->STRUCTURES (this finds all the candidate structures for a gene and the coverage)
            (currently it also makes a pdb input as per the last stage but it will be redone when stage 3 is completed)
                                |                       
            *SCRIPT 03: STRUCTURES->STRUCTURE (NOT DONE) this aggregates the structures into 1/makes a choice
                                |                       
            *SCRIPT 04: FOLDX PREPARE (This prepares the inputs for a foldx batch)
                                                           
 ```
-----------------------------------------------------------------------
### Overview
This pipeline manipulates data from above and ends up with some protein structures for a gene under selection, and for each one there is a list of variants which can be fed into the foldx pipeline.

-----------------------------------------------------------------------
### Inputs

1. We do not yet have the first script, but the inputs are probably a vcf file
2. The vcf file uses an r script dNdS to create 2 files, a list of variants and a list of genes under selection, these feed into 02
3. The list of structures and coverage is used to make a final structure decision
4. The final structure is used to create the inputs for the foldx pipeline

-----------------------------------------------------------------------
### Pipeline parent script

- Only the name of the dataset is needed to run the pipeline, e.g. we have the name "shearwater" for a test set
-----------------------------------------------------------------------
### Scripts
- pipeline02_genestoproteins: takes the lists from the dataset dNdS output and turns it into candidate pdb files with variant coverage


-----------------------------------------------------------------------
### Structure
For each py there is an sh file so the batch can be run from either python or command line-hpc

-----------------------------------------------------------------------





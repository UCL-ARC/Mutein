##### Owner:Rachel
-----------------------------------------------------------------------
## HPC FoldX Pipeline for myriad
-----------------------------------------------------------------------
```
The scripts structure is:

                                    SCRIPT 01:FOLDX REPAIR
                                        |                |
                            SCRIPT 02:SPLIT          SCRIPT 05:VARIANT SPLIT
                                        |                |
            [ARRAY JOBS]    SCRIPT 03:FOLDX POSSCAN  SCRIPT 06:FOLDX BUILD
                                        |                |
                            SCRIPT 04:AGGREGATE      SCRIPT 07:VARIANT AGGREGATE
 ```
-----------------------------------------------------------------------
### Overview
This pipeline compares the background protein folding stability of a given pdb structure against known variants.

- The background folding stability is calculated by mutating every single possible residue in a protein structure
- The variants are calculating by perumtating the possibilities of the known variants
-----------------------------------------------------------------------
### Inputs
1. A pdb file
2. A variants file
-----------------------------------------------------------------------
### Pipeline parent script
- The 00 script is a parent script that runs all the other scripts with the correct dependencies
-----------------------------------------------------------------------
### Scripts
1. This runs foldx repairs on a pdb file, the repairs ensure that the atom positions are relaxed into a favoured position.
- background mutations
2. This splits the possible mutations up into a config file ready for parallising
3. An array job: This calls foldx positionscan function with 1 of the rows from the config file, 1 set of mutations
4. This aggregates and perfoms some analysis on all the individual muations
- variant mutations
5. This permutes the possible mutations and makes a config file
6. An array job: This calls build for each mutation
7. And finally aggregates and analyses

- N.B. there is not yet a job that compares the background to the variant data
- N.B.2 the number of array jobs is arbitrary for script 3 and specific to the number of variants for script 6
-----------------------------------------------------------------------
### Structure
There are pairs of scripts bash/python. The bash submits to qsub, the python runs in python. This enables a test environment to run the python scripts via both CI and user testing from the tests directory.

There are 3 levels of paramaters. In order of least preferred to most:
- Batch params (eg number of array jobs, time allowed, script dependence)
- Pdb parameters (eg override of arrays for variants, name of variant, number of repairs)
- Command line inputs - pdb needed, override anything else if you want (eg change the number of array jobs for testing)
-----------------------------------------------------------------------
### Environment
The scripts can be run in 3 modes - python, hpc and inputs (inputs only records the inputs to external libraries). You can override this with user= in the batch.
You can add your dev invironment to the top of helper.py
An example of a batch command with parameters overridden:



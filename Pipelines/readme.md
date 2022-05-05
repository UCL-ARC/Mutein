
# Pipeline Overview
#### 1)
#### 2)
#### 3) GeneProt: This starts with I think a vcf file and ends with the protein structures
#### 4) FoldX: This pipeline takes a pdb file and a list of mutations and creates data frames of ddg
#### 4) GeneStitch: This pipeline puts all the structures back together to show ddg for a gene



# Instructions for users


# Instructions for developers
There is a framework set up so it should be relatively easy to create new pipelines that can be submitted to qsub or run locally.
- Each pipeline needs a yaml script to set it up, and a batch sh script to run on the server.
- Additionally the way the batch runs expects the python script parameters tto be in a certain format.
- Best practice: also add a test to the tests directory.
More information is in the detailed document: [Creating a new pipeline](../Documents/new_pipeline.md)




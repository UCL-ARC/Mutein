#### Owner:Rachel

### HPC FoldX Pipeline for myriad

#### Pipeline inputs
##### A pdb file
##### A list of mutations in a file of the form
Mutation,Variant
NA501Y,Alpha
AA570D,Alpha
TA716I,Alpha
SA982A,Alpha
DA1118H,Alpha
DA614G,Alpha

#### The data is structured in folders
inputs/setname/*
scripts/*
results/setname/*
outputs/setname*
(The results and outputs are not are not synced with github)

To run the pipeline on myriad:
clone the repo
#### >> module load python3/recommended
#### >> python3 scripts/Sh00_Myriad_pipeline.py jobs=1234567 rows=45 setname=6vxx_45 length=6:00:0 (or whatever parameters)

#### All python scripts take the same consistent parameters, and can be run via this hpc batch mechanism, or locally:


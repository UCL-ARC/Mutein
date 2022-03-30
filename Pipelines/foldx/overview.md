##### Owner:Rachel
########################################################################
## HPC FoldX Pipeline for myriad
########################################################################
### 1.Pipeline inputs
##### A pdb file
##### A list of mutations in a file of the form
Mutation,Variant
NA501Y,Alpha
AA570D,Alpha
TA716I,Alpha
SA982A,Alpha
DA1118H,Alpha
DA614G,Alpha
########################################################################
### 2.The data is structured in folders
scripts/*
inputs/pdb/*
thruputs/pdb/*
results/setname/*
outputs/setname/*
(The thruputs, results and outputs are not are not synced with github)
########################################################################
### 3.To run the pipeline on myriad:
#### clone the repo
########################################################################
#### >> module load python3/recommended
#### >> python3 scripts/Sh00_Myriad_pipeline.py jobs=1234567 rows=45 setname=6vxx_45 length=6:00:0 (or whatever parameters)
########################################################################
#### The batch can be run manually and on each individial script too, environment variables 
#### (you may need to add yoruself) will recognise whether you are running locally or on a server
########################################################################


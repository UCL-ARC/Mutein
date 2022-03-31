##### Owner:Rachel
########################################################################
## HPC FoldX Pipeline for myriad
########################################################################
### 0. When everthing is installed and running, you just need to do the following to run the pipeline
```
>> cd ~/MuteinPipeline/foldx/scripts/
>> module load python3/recommended
>> python3 foldx00_pipeline.py jobs=1234567 split=45 setname=6vxx_45 length=6:00:0 (or whatever parameters)
```
### TODO ###
1. It has my ucl username embedded in the bash scripts for qsub working directory
2. I force everyone to use the same path to install (well, you can change the install script)

#### -- But to get to that state - install or update as follows:
### 1.Pipeline inputs
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
########################################################################
### 2.The data is structured in folders
```
scripts/*
inputs/pdb/*
thruputs/pdb/*
results/setname/*
outputs/setname/*
(The thruputs, results and outputs are not are not synced with github)
```
########################################################################
### 3.To run the pipeline on myriad:
##### clone the repo
```
# 1. Get the code
>> git checkout .
>> git pull
>> cd Mutein
# 2. Run the install script, let's not run everything within github on the servers(it is ok on your local machine)
>> chmod +x Pipelines/foldx/install.sh
>> Pipelines/foldx/install.sh
# 3. This has created/copied the scripts in ~/MuteinPipeline/foldx/scripts/
# Now navigate there and run the pipeline
>> cd ~/MuteinPipeline/foldx/scripts/
>> module load python3/recommended
>> python3 foldx00_pipeline.py jobs=1234567 split=45 setname=6vxx_45 length=6:00:0 (or whatever parameters)
```
########################################################################
##### The batch can be run manually and on each individial script too, environment variables will recognise
##### whether you are running locally or on a server (you may need to add yourself)
########################################################################

### Continuous integration



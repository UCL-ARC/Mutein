##### Owner:Rachel
--------------------------------------------------------------------------
## HPC FoldX Pipeline for myriad
--------------------------------------------------------------------------
### When everthing is installed and running, you just need to do the following to run the pipeline
```
>> cd ~/MuteinPipeline/foldx/scripts/
>> module load python3/recommended
>> python3 foldx00_pipeline.py jobs=1234567 pdb=6vxx
```
--------------------------------------------------------------------------
#### -- But to get to that state - 
--------------------------------------------------------------------------
### 1.The data is structured in folders
```
scripts/*
inputs/pdb/*
thruputs/pdb/*
results/setname/*
outputs/setname/*
(The thruputs, results and outputs are not are not synced with github)
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
# Run the install script
>> chmod +x Pipelines/foldx/install.sh
>> Pipelines/foldx/install.sh
# 3. This has created/copied the scripts in ~/MuteinPipeline/foldx/scripts/
# Now navigate there and run the pipeline
>> cd ~/MuteinPipeline/foldx/scripts/
>> module load python3/recommended
>> python3 foldx00_pipeline.py jobs=1234567 pdb=6vxx (or whatever parameters)

```
-----------------------------------------------------------------------
```
The scripts dependency is:

                                    SCRIPT 01: FOLDX REPAIR
                                        |                |
                            SCRIPT 02: SPLIT          SCRIPT 05: VARIANT SPLIT
                                        |                |
            [ARRAY JOBS]    SCRIPT 03: FOLDX POSSCAN  SCRIPT 06: FOLDX BUILD
                                        |                |
                            SCRIPT 04: AGGREGATE      SCRIPT 07: VARIANT AGGREGATE
 ```
-----------------------------------------------------------------------
##### The batch can be run manually and on each individial script too, environment variables will recognise
##### whether you are running locally or on a server (you may need to add yourself)
-------------------------------------------------------------------------------------------------
### Useful docs
### Continuous integration
https://www.techiediaries.com/python-unit-tests-github-actions/

##### Using VSCode and WSL
- necessarily easy (access and security problems, but the beginning of this works)Not
-- Virtual env in vscode: https://techinscribed.com/python-virtual-environment-in-vscode/#:~:text=Using%20Python%20Virtual%20Environment%20in%20VSCode%201%20Install,installed%2C%20VSCode%20will%20show%20an%20error%20like%20this.
- In vscode once created, the venv can be selected in the bottom-right next to "Python3"
- Install things into the venv from the command line within VSCode 






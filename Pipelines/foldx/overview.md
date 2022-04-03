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
# Clone the main branch or a branch
>> cd Mutein

# 1. Get the code
>> git checkout .
>> git pull
# 2. Run the install script, let's not run everything within github on the servers(it is ok on your local machine)
>> chmod +x Pipelines/foldx/install.sh
>> Pipelines/foldx/install.sh
# 3. This has created/copied the scripts in ~/MuteinPipeline/foldx/scripts/
# Now navigate there and run the pipeline
>> cd ~/MuteinPipeline/foldx/scripts/
>> module load python3/recommended
>> python3 foldx00_pipeline.py jobs=1234567 split=45 setname=6vxx_45 length=6:00:0 (or whatever parameters)

Stabdard covid call is
python3 foldx00_pipeline.py jobs=1234 id=3@array=50 id=3@time=6:00:0 setname=6vxx50 repairs=5
python3 foldx00_pipeline.py jobs=1234 id=3@array=200 id=3@time=3:00:0 setname=6vxx200 repairs=10

```
########################################################################
##### The batch can be run manually and on each individial script too, environment variables will recognise
##### whether you are running locally or on a server (you may need to add yourself)
########################################################################

### Continuous integration



### Useful docs
1. Virtual env in vscode: https://techinscribed.com/python-virtual-environment-in-vscode/#:~:text=Using%20Python%20Virtual%20Environment%20in%20VSCode%201%20Install,installed%2C%20VSCode%20will%20show%20an%20error%20like%20this.




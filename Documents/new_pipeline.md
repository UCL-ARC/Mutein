
#### (RSA 5.5.22)

# Mutein Documents: Creating a New Pipeline

There is a framework set up so it should be relatively easy to create new pipelines that can be submitted to qsub or run locally.

## Python script
The main work is in a python script (R scripts to be added later)
The python scripts take a specific format in order to process the inputs from the batch system:
```
    def run_pipeline01(args):
```
where args[1] is a long string of inputs in the format key=value@key=value@...e.g.
```
    pdb=6vxx@repairs=5
```
The arguments can be parsed out easily with the arguments class:
```
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
```
## YAML configuration of the batch
The batch is configured with a yaml file which poinys to the python files above. They are executed in order,a nd where run through qsub, dependencies can be set. The batch can be run on python mode too for debugging.
Example batch config:
```
---
id: 'Foldx1'
qsub_id: 'foldx_repair'
pipe_dir: 'foldx/scripts/'
script: 'foldx01_repair'
time: '3:00:0'
dependency: '-1'
array: 0
inputs: "repairs=5"
active: 'Y'
---
id: 'Foldx2'
qsub_id: 'foldx_params'
work_dir: 'foldx/scripts/'
script: 'foldx02_makeparams'
time: '1:00:0'
dependency: 'Foldx1'
array: 0
inputs: "repairs=5@split=100"
active: 'Y'
---
```
- Use the id for dependencies, other qsub specific params are the time and memory. 
- work_dir is the directory where the script lives.
- array is for array jobs, if it is other than 0 and array batch is submitted
- active set to N means it is skipped(testing)
- qsub_id is what you will see in the qsub admin

## Array jobs
When you have an array job, all that changes is that the array job will pass into the python script the task number appended to the unputs, so for example the inputs above would become:
```
    pdb=6vxx@repairs=5@task=1
```
You can handle the task number however you like in the python script.
The template bash script for array jobs and single jobs are:
- pipeline_array.sh
- pipeline_single.sh
WARNING these should not be changed (unless fixed) as everything uses them, they are templates. Change your python scripts.

## Bash script for QSUB
The bash script is simply a call to the template python script that knows how to submit to qsub.
The template script is called pipeline_qsubber.py and doesn't need changing.... unless....

...The bash scripts match up the template script and yaml files along with 3 specific parameters for the domain:

python pipeline_qsubber.py batch_pdb.yml qsub "dataset" "gene" "pdb" "chain"
e.g.```
python pipeline_qsubber.py batch_pdb.yml qsub notch NOTCH1 AF-P46531-F1-model_v2 A
```
The reason for these is that they provide a template to the yaml files so that a single yaml file can be used for different structures.
The inputs are appended to the "inputs" with the same format pdb= @ etc
**TODO so I see now that I could have just passed params in that format in and it would be entirely flexible**

## Testing
Each python script can be automatically run for debugging as per below.
```
retpath = "/".join(dirs) + "/genestitch"
sys.path.append(retpath)
import Paths

def addpath():        
    print(sys.path)

def test_geneprot_pipeline(inputs):
    addpath()
    import pipeline02_genestoproteins as p02
    args = ["", inputs]
    p02.run_pipeline(args)

test_geneprot_pipeline("dataset=notch@")
```


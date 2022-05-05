
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
    pdb=6vxx@repairs=5@task=1
```
The arguments can be parsed out easily with the arguments class:
```
    argus = Arguments.Arguments(args)
    dataset = argus.arg("dataset")
```


- Each pipeline needs a yaml script to set it up, and a batch sh script to run on the server.
- Additionally the way the batch runs expects the python script parameters tto be in a certain format.
- Best practice: also add a test to the tests director


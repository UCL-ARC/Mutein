"""

------------------------
RSA 29/03/22
------------------------
Pipeline overview script for Mutein
       PIPELINE01 - GENES TO PROTEINS
       PIPELINE02 - FOLDX REPAIR
                                    
------------------------

This pipline script is called with a single input of a batch config yaml file
It takes 2 inputs
- the batch.yml file
- qsub, sh or py: where qsub means it runs the files on qsub and py means it runs the files rough python directly and sh the sh files directly
The yaml contains the scripts and dependencies
This is supposed to be a pretty straightforward file that doesn't need to be looked at much, the config defines the batch
---
id: '1'
script: 'geneprot/scripts/pipeline02_genestoproteins'
time: '3:00:0'
dependency: '-1'
array: 0
inputs: dataset=shearwater
---
"""
import os
import pwd
import yaml

#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Config
import Paths
import Arguments
import QSubRunner as qsub
import SubRunner as sub


##### INPUTS #############################################
## Pipeline jobs sequence

def overall_rsa(args):
    # The environment is loaded by arguments by default
    argus = Arguments.Arguments([])
    ret_array = []
    print("#### FOLDX PIPELINE - batch creation ####")
    # There are only 2 arguments
    yaml_file = args[1] #1) a yaml file path with the batch definition    
    py_or_sh = args[2] #2) qsub or py or sh for python or hpc batch or just sh
    # We want the user
    homeuser = pwd.getpwuid(os.getuid())[0]
    print("HomeUser=", homeuser)
    print("Running",yaml_file," for user=",homeuser)
    # We want to be in the script directory    
    script_path = os.path.dirname(os.path.realpath(__file__))
    dir_path =  script_path + "/"
    print("# overall pipeline: changing directory to", dir_path)
    os.chdir(dir_path)

    # Load in the batch yaml file, which contains ONLY batch related parameters
    # apart from the final line inputs which are space delim inputs for the script
    batch_dic = {}
    batch_list = []
    with open(yaml_file, "r") as fr:
        pipes = yaml.safe_load_all(fr)
        for pipe in pipes:
            if pipe != None:
                print('pipe|',pipe)
                id = pipe["id"]
                script = pipe["script"]
                time = pipe["time"]
                dependency = pipe["dependency"]
                array = pipe["array"]
                inputs = pipe["inputs"]
                batch_dic[str(id)] = (script, time, dependency, array, inputs)
                batch_list.append(str(id))

    dependencies = {}
    for id in batch_list:        
        script, time, dependency, array, inputs = batch_dic[id]
        print("# overall pipeline script:",id,script, time, dependency, array, inputs)        
        if py_or_sh == "qsub":
            if dependency != -1:
                if dependency in dependencies:                
                    dependency = dependencies[id]
                else:
                    dependency = -1
            runner = qsub.QSubRunner(script + ".sh", dependency, time, array,homeuser,inputs)
            dep = runner.run()
            dependencies[id] = dep
        elif py_or_sh == "py":
            runner = sub.SubRunner(argus.arg("pythonexe"), script + ".py", inputs)
            dep = runner.run()
        elif py_or_sh == "sh":                                    
            dirs = (script_path + "/" + script).split("/")[:-1]
            newpath = "/".join(dirs)                     
            print("# overall pipeline: changing directory to", newpath)
            os.chdir(newpath)
            runner = sub.SubRunner("", script_path + "/" + script + ".sh", inputs)
            dep = runner.run()


    
    

####################################################################################################
if __name__ == "__main__":
    import sys
    globals()["overall_rsa"](sys.argv)


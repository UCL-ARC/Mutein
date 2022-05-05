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

def pipeline_qsubber(args):
    # The environment is loaded by arguments by default
    print("ARGS=",args)
    argus = Arguments.Arguments(args,spaced=False)
    

    # There are 6 arguments
    yaml_file = args[1] #1) a yaml file path with the batch definition    
    py_or_sh = args[2] #2) qsub or py or sh for python or hpc batch or just sh
    # everything is defined in the yaml APART from dataset, gene, pdb        
    dataset,gene,pdb,chain = "","","",""
    if len(args) > 2:
        dataset = args[3] #3) dataset
    if len(args) > 3:
        gene = args[4] #4) gene
    if len(args) > 4:
        pdb = args[5] #5) pdb
    if len(args) > 5:
        chain = args[6] #6) chain
    print("#### MUTEIN PIPELINE ####", yaml_file,py_or_sh, dataset, gene, pdb,chain)

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
                print('### yaml load|',pipe)
                id = pipe["id"]
                qsubid = pipe["qsub_id"]
                work_dir = pipe["work_dir"].strip()
                script = pipe["script"].strip()
                time = pipe["time"].strip()
                dependency = pipe["dependency"].strip()
                array = pipe["array"]
                inputs = pipe["inputs"].strip()
                active = pipe["active"].strip()=="Y"
                # add the dataset, gene and pdb onto inputs
                if len(inputs) > 0:
                    inputs += "@"
                inputs += "dataset=" + dataset
                inputs += "@gene=" + gene
                inputs += "@pdb=" + pdb
                inputs += "@chain=" + chain
                if active:
                    batch_dic[str(id)] = (qsubid,work_dir,script, time, dependency, array, inputs)
                    batch_list.append(str(id))

    dependencies = {}
    for id in batch_list:
        isarray =  int(array)>0
        qsubid,work_dir,script, time, dependency, array, inputs = batch_dic[id]
        #print("# overall pipeline script:",id,qsubid,work_dir,script, time, dependency, array, inputs)        
        if "qsub" in py_or_sh:
            dep = -1
            if str(dependency) != "-1":
                if dependency in dependencies:                
                    dep = dependencies[dependency]
                else:
                    dep = -1            
            runner = qsub.QSubRunner(id,qsubid,script, dir_path,work_dir,dep,time,array,homeuser,inputs,py_or_sh!="qsub")
            dep = runner.run()
            dependencies[id] = dep
        elif py_or_sh == "py":
            runner = sub.SubRunner(argus.arg("pythonexe"),dir_path,work_dir,script,".py",inputs,isarray)
            dep = runner.run()
        elif py_or_sh == "sh":                                                
            runner = sub.SubRunner("bash",dir_path,work_dir,script,".sh",inputs,isarray)
            dep = runner.run()


    
    

####################################################################################################
if __name__ == "__main__":
    import sys
    globals()["pipeline_qsubber"](sys.argv)


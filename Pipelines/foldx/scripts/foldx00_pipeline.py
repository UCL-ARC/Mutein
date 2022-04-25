"""
------------------------
RSA 29/03/22
------------------------
Pipeline overview script for foldx job on Myriad

This script is a parent script that runs the entire foldx pipeline at a top level.
It manages the dependencies between the scripts and allows you to run as either HPC, python or inputs mode

The scripts dependency is:

                                    SCRIPT 01: FOLDX REPAIR
                                    |                |
                        SCRIPT 02: SPLIT     SCRIPT 05: VARIANT SPLIT
                                |                       |
        [ARRAY JOBS]    SCRIPT 03: FOLDX POSSCAN     SCRIPT 06: FOLDX BUILD
                                |                       |
                        SCRIPT 04: AGGREGATE        SCRIPT 07: VARIANT AGGREGATE
------------------------
"""
import os
import pwd
import subprocess
import pandas as pd
import yaml

#import from the shared library in Mutein/Pipelines/shared/lib
import sys
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + '/shared/libs'
sys.path.append(retpath)
import Config
import Paths
import Arguments


##### INPUTS #############################################
## Pipeline jobs sequence
# 1= repairing pdb
# 2= making param file for splits
# 3= performing the position scan ddg mutations (parallel)
# 4= aggregating 3
# 5= making params file for variants
# 6 = performing variant ddg (parallel)
# 7 = aggregating 6


def run_pipeline00(args):
    ret_array = []
    print("#### FOLDX PIPELINE - batch creation ####")
    #The arguments HAVE to include a pdb name
    argus = Arguments.Arguments(args)    
    pdbcode = argus.arg("pdb")
    pdb_path = Paths.Paths("pdb",dataset="",gene="",pdb=pdbcode)    
    pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    argus.addConfig(pdb_config.params)

    ### Change into script directory
    dir_path = os.path.dirname(os.path.realpath(__file__)) + "/"
    print("## ... changing directory to", dir_path)
    os.chdir(dir_path)
    homeuser = pwd.getpwuid(os.getuid())[0]
    print("HomeUser=", homeuser)
    ### Process paramaters in order of preference, job, config, pipeline
    
    #cfgplparams = hlp.configpipelineparams(argus.arg("pdb"))
    pipelineparams = argus.pipelineparams
    #print("Pipelines=", pipelineparams)
    # script extension is either sh or py depending on bash or python environment
    ext = ".sh"
    if argus.arg("env") == "python":
        ext = ".py"
    # The batch is defined in the file batch.yml
    batch_dic = {}
    with open("batch.yml", "r") as fr:
        pipes = yaml.safe_load_all(fr)
        for pipe in pipes:
            # print('pipe|',pipe)
            id = pipe["id"]
            script = pipe["script"]
            time = pipe["time"]
            dependency = pipe["dependency"]
            array = pipe["array"]
            batch_dic[str(id)] = (script, time, dependency, array)
    #############################################################
    dependencies = {}
    runs = []
    for j in argus.arg("jobs"):
        if str(j) in batch_dic:
            script, time, dependency, array = batch_dic[str(j)]
            # check there are no overrides from inuts
            if j in pipelineparams:
                idparams = pipelineparams[j]
                if "time" in idparams:
                    time = idparams["time"]
                if "array" in idparams:
                    array = idparams["array"]
            dep = "-1"
            if str(dependency) != "-1" and str(dependency) in argus.arg("jobs"):
                dep = dependency
            runs.append([j, "qsub", script + ext, dep, time, array])
            dependencies[str(j)] = str(dep)

    for job, exe, script, dependency, time, array in runs:
        # print(job,exe,script,dependency)
        if argus.arg("env") == "hpc":
            os.system("chmod +x " + script)
        args = []
        if "python" in argus.arg("env"):
            args.append(argus.arg("pythonexe"))
        else:
            args.append("qsub")
            if str(dependency) != "-1":
                args.append("-hold_jid")
                args.append(dependency)
            if int(array) > 0:
                args.append("-t")
                args.append("1-" + str(array))
            args.append("-l")  # $ -l h_rt=5:00:0
            args.append("h_rt=" + time)
            args.append("-wd")  # $ -wd /home/ucbtlcr/Scratch/workspace
            args.append(
                "/home/" + homeuser + "/Scratch/workspace"
            )  # $ -wd /home/ucbtlcr/Scratch/workspace

        args.append(script)
        if (
            argus.arg("env") == "hpc"
            or argus.arg("env") == "inputs_hpc"
        ):
            args.append(argus.arg("pdb"))  # 1
            args.append("xxx")  # 2
            args.append(argus.arg("split"))  # 3
            args.append(argus.arg("mutation"))  # 4
            args.append(argus.arg("variant"))  # 5
            args.append("variants.csv")  # 6
            args.append(argus.arg("repairs"))  # 7
        else:
            args.append(script)
            args.append("pdb=" + argus.arg("pdb"))  # 1            
            args.append("name=xxx")  # 2
            args.append("split=" + str(argus.arg("split")))  # 3
            args.append("mutation=" + argus.arg("mutation"))  # 4
            args.append("variant=" + argus.arg("variant"))  # 5
            args.append("variantfile=variants.csv")  # 6
            args.append("repairs=" + str(argus.arg("repairs")))  # 7

        print(args)
        ret_array.append(args)
        if argus.arg("env") == "hpc":
            # print('Running on hpc')
            process = subprocess.Popen(
                args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            result = process.communicate()
            print(result)  # e.g. Your job 588483 ("foldx-aggregate") has been submitted
            results = result[0].split(" ")
            jobid = results[2]
            if "." in jobid:
                results = jobid.split(".")
                jobid = results[0]
            dependencies[int(job)] = jobid
        elif argus.arg("env") == "python":
            # print('Running in python')
            dependencies[int(job)] = 0
            process = subprocess.Popen(
                args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            result = process.communicate()
            print(result)  # e.g. Your job 588483 ("foldx-aggregate") has been submitted
        else:
            # print('Not running')
            dependencies[int(job)] = 0

    return ret_array


####################################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline00"](sys.argv)

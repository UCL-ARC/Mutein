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
from os.path import exists

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")
retpath = "/".join(dirs) + ""
sys.path.append(retpath)
import Arguments
import QSubRunner as qsub
import SubRunner as sub
import FileDf
import Paths

##### INPUTS #############################################
## Pipeline jobs sequence


def pipeline_qsubber(args):
    # The environment is loaded by arguments by default
    print("ARGS=", args)
    argus = Arguments.Arguments(args, spaced=False)

    ## $PWD ${config} $run notch NOTCH1 AF-P46531-F1-model_v2

    # There are 7 arguments
    install_dir = args[
        1
    ]  # 1) the executable installation directory, the root directory of the peipeline
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    working_dir = args[
        2
    ]  # 1) the working dir, the root that the data output and input lives in
    yaml_file = args[3]  # 2) a yaml file path with the batch definition
    py_or_sh = args[4]  # 3) qsub or py or sh for python or hpc batch or just sh
    # everything is defined in the yaml APART from 3 template inputs
    dataset, gene, pdb = "", "", ""
    dataset = args[5]  # 3) dataset
    gene = args[6]  # 4) gene
    pdb = args[7]  # 4) pdb
    if dataset == "x":
        dataset = ""
    if gene == "x":
        gene = ""
    if pdb == "x":
        pdb = ""

    # last_stat = args[8]  # 5) the last status or count
    # if last_stat == "x":
    #    print("!!! qsub in fail state, exiting")
    #    return "x"

    # count = int(last_stat)
    count = 0

    print(
        "#### MUTEIN PIPELINE ####",
        install_dir,
        working_dir,
        yaml_file,
        py_or_sh,
        dataset,
        gene,
    )

    names_and_ids = []

    # There are 3 files that can control the number of tasks to be run. If those files exists we load up those numbers now.
    gene_tasks = []
    pdb_tasks = 0
    pdb_tasks_missing = 0
    params_tasks = 0
    vparams_tasks = 0
    unparams_tasks = 0
    vunparams_tasks = 0

    if gene == "ALL":
        gene = ""
        path = Paths.Paths(
            working_dir, install_dir, dataset=dataset
        )
        gene_tasks_file = path.inputs + "genes_pdb_list.csv"
        print("Gene file=", gene_tasks_file)
        if exists(gene_tasks_file):
            with open(gene_tasks_file) as fr:
                lines = fr.readlines()
                for l in range(1, len(lines)):
                    line = lines[l]
                    if "," in line:
                        gene = line.split(",")[1].strip()
                        gene_tasks.append(gene)
        else:
            gene_tasks = [""]
    else:
        gene_tasks = [gene]

    for gene in gene_tasks:
        print(
            "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Creating for gene=",
            gene,
            "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",
        )
        path = Paths.Paths(
            working_dir,
            install_dir + "Pipelines/geneanalysis",
            dataset=dataset,
            gene=gene,
            pdb=pdb,
        )
        pdb_tasks_file = path.outputs + "pdb_tasklist.csv"
        if exists(pdb_tasks_file):
            with open(pdb_tasks_file) as fr:
                lines = fr.readlines()
                pdb_tasks = len(lines) - 1
        pdb_tasks_file_missing = path.outputs + "pdb_tasklist_incomplete.csv"
        if exists(pdb_tasks_file_missing):
            with open(pdb_tasks_file_missing) as fr:
                lines = fr.readlines()
                pdb_tasks_missing = len(lines) - 1
        params_tasks_file = path.inputs + "params_background.txt"
        if exists(params_tasks_file):
            with open(params_tasks_file) as fr:
                lines = fr.readlines()
                params_tasks = len(lines) - 1
        unparams_tasks_file = path.inputs + "params_background_incomplete.txt"
        if exists(unparams_tasks_file):
            with open(unparams_tasks_file) as fr:
                lines = fr.readlines()
                unparams_tasks = len(lines) - 1
        vparams_tasks_file = path.inputs + "params_variants.txt"
        if exists(vparams_tasks_file):
            with open(vparams_tasks_file) as fr:
                lines = fr.readlines()
                vparams_tasks = len(lines) - 1
        vunparams_tasks_file = path.inputs + "params_variants_incomplete.txt"
        if exists(vunparams_tasks_file):
            with open(vunparams_tasks_file) as fr:
                lines = fr.readlines()
                vunparams_tasks = len(lines) - 1

        # We want the user
        homeuser = pwd.getpwuid(os.getuid())[0]
        print("HomeUser=", homeuser)
        print("Running", yaml_file, " for user=", homeuser)
        # We want to be in the script directory
        # script_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = install_dir + "Pipelines/"
        # print("# overall pipeline: changing directory to", dir_path)
        # os.chdir(dir_path)

        # Load in the batch yaml file, which contains ONLY batch related parameters
        # apart from the final line inputs which are space delim inputs for the script
        batch_dic = {}
        batch_list = []
        with open(yaml_file, "r") as fr:
            pipes = yaml.safe_load_all(fr)
            for pipe in pipes:
                if pipe != None:
                    print("### yaml load|", pipe)
                    id = pipe["id"]
                    qsubid = pipe["qsub_id"]
                    pipe_dir = pipe["pipe_dir"].strip()
                    script = pipe["script"].strip()
                    time = pipe["time"].strip()
                    dependency = pipe["dependency"].strip()
                    cores = pipe["cores"].strip()
                    array = pipe["array"]
                    if int(array) == -1:
                        arrayfile = pipe["arrayfile"].strip()
                        if arrayfile == "pdbs":
                            array = pdb_tasks
                        if arrayfile == "unpdbs":
                            array = pdb_tasks_missing                            
                        if arrayfile == "params":
                            array = params_tasks
                        if arrayfile == "vparams":
                            array = vparams_tasks
                        if arrayfile == "unparams":
                            array = unparams_tasks
                        if arrayfile == "vunparams":
                            array = vunparams_tasks
                        if int(array) == 0:
                                print("Complete",gene)
                                continue
                        
                    inputs = pipe["inputs"].strip()
                    active = pipe["active"].strip() == "Y"
                    # add the dataset, gene and pdb onto inputs
                    if len(inputs) > 0:
                        inputs += "@"
                    inputs += "install_dir=" + install_dir
                    inputs += "@data_dir=" + working_dir
                    inputs += "@dataset=" + dataset
                    inputs += "@gene=" + gene
                    inputs += "@pdb=" + pdb
                    id_ext = ""
                    if len(dataset)>=2:
                        id_ext += dataset[0:2]
                    if len(gene)>=2:
                        id_ext += ""+gene[0:2]
                    if len(pdb)>=2:
                        id_ext += ""+pdb[0:2] + "_"
                    qsubid = id_ext + qsubid

                    # inputs += "@chain=" + chain
                    if active:
                        batch_dic[str(id)] = (
                            qsubid,
                            pipe_dir,
                            script,
                            time,
                            dependency,
                            array,
                            inputs,
                            cores,
                        )
                        batch_list.append(str(id))

        dependencies = {}
        for id in batch_list:
            (
                qsubid,
                pipe_dir,
                script,
                time,
                dependency,
                array,
                inputs,
                cores,
            ) = batch_dic[id]
            isarray = int(array) > 0
            # print("# overall pipeline script:",id,qsubid,pipe_dir,script, time, dependency, array, inputs)
            if "qsub" in py_or_sh:
                dep = -1
                if str(dependency) != "-1":
                    if dependency in dependencies:
                        dep = dependencies[dependency]
                    else:
                        dep = -1
                runner = qsub.QSubRunner(
                    id,
                    qsubid,
                    script,
                    install_dir,
                    working_dir,
                    pipe_dir,
                    dep,
                    time,
                    array,
                    homeuser,
                    inputs,
                    cores,
                    py_or_sh != "qsub",
                )
                dep = runner.run()
                if dep == "x":
                    print("!!!Abandoning submissions!!!")
                    return "x"
                else:
                    count += 1
                    dependencies[id] = dep
                    names_and_ids.append([qsubid, dep])
            elif py_or_sh == "py":
                runner = sub.SubRunner(
                    argus.arg("pythonexe"),
                    install_dir,
                    working_dir,
                    pipe_dir,
                    script,
                    ".py",
                    inputs,
                    isarray,
                )
                dep = runner.run()
            elif (
                py_or_sh == "sh"
            ):  # TODO make this go simply straight through passing all inputs
                runner = sub.SubRunner(
                    "bash",
                    install_dir,
                    working_dir,
                    pipe_dir,
                    script,
                    ".sh",
                    inputs,
                    isarray,
                )
                dep = runner.run()

    return str(count)


####################################################################################################
if __name__ == "__main__":
    import sys

    globals()["pipeline_qsubber"](sys.argv)

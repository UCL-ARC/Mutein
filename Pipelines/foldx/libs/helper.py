"""
------------------------
RSA 29/03/22
------------------------
Helper module for pipeline script for foldx job on Myriad
----
"""

import os

### These functions consistently handle the paramater inputs for the script, merging config and overrides
def addlinetoparams(arg, params):
    args = []
    if "=" in arg and "@" not in arg:
        args = arg.split("=")
        p, v = args[0], args[1]
        params[p] = v
    return params


def addpipelinetoparams(arg, params):
    args = []
    if "=" in arg and "@" in arg:
        args = arg.split("@")
        print(args)
        id = args[0].split("=")[1]
        pv = args[1].split("=")
        p, v = str(pv[0]), str(pv[1])
        print(args[0], args[1], id, p, v)
        if id not in params:
            params[id] = {}
        params[id][p] = str(v)
    return params


def configparams(pdb):
    # set up some defaults for any batch to run without paramaters
    params = {}
    params["jobs"] = "1234567"
    params["chain"] = "A"
    params["pdb"] = "6vxx"
    params["name"] = "6vxx_50"
    params["row"] = "1"
    params["mutation"] = "."
    params["time"] = "."
    params["variant"] = "Alpha"
    if pdb != "":
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = dir_path[:-7]
        input_path = dir_path + "inputs/"
        configfile = input_path + pdb + "/config.cfg"
        with open(configfile) as fr:
            cfgcontent = fr.readlines()
            for line in cfgcontent:
                line = line.strip()
                params = addlinetoparams(line, params)
    if "variantfile" not in params:
        params["variantfile"] = params["pdb"] + "_vars"

    # Code to decide if it is test or live environment can be overridden from the command line
    import environments

    envs = environments.getenvironment()
    params["user"] = envs["user"]
    return params


def configpipelineparams(pdb):
    # set up some defaults for any batch to run without paramaters
    params = {}
    if pdb != "":
        dir_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = dir_path[:-7]
        input_path = dir_path + "inputs/"
        configfile = input_path + pdb + "/config.cfg"
        with open(configfile) as fr:
            cfgcontent = fr.readlines()
            for line in cfgcontent:
                line = line.strip()
                params = addpipelinetoparams(line, params)
    return params


def pipelineparams(argvs, params):
    for i in range(1, len(argvs)):
        arg = argvs[i]
        params = addpipelinetoparams(arg, params)
    return params


def mergeparams(configparams, jobparams):
    # the job params take precendence
    for cfg, val in jobparams.items():
        configparams[cfg] = val
    return configparams


def goto_job_dir(dir_path, args, params, name):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    os.chdir(dir_path)
    inputs_file = name + ".log"
    with open(inputs_file, "w") as fw:
        for arg in args:
            fw.write(str(arg) + " ")
        fw.write("\n")
        for cfg, val in params.items():
            fw.write(str(cfg) + "=" + str(val) + "\n")

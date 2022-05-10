"""
RSA 5/4/22
------------------------
Class to manage the arguments and precendence of the argument layer
The command line overrides everything
The config file overrides the batch
The batch params ar ethe first level of configutation (script order, dependecies, time, array splits)

This class is an effort to standardise the inputs across pipelines
The first version of this is in the foldx pipeline, it will need to be migrated to this which
as I have gne up the pipeline reflects a new more generic structure (RA)
"""

import environments


class Arguments:
    def __init__(self, inputs, spaced=True):
        # combine config and job input params
        print("### Arguments", inputs)
        if spaced:
            self.params = self.spacedparams(inputs)
        else:
            self.params = self.inputparams(inputs)
            self.pipelineparams = self.inputparams(inputs)
        if "user" in self.params:
            envs = environments.getenvironment(self.params["user"])
        else:
            envs = environments.getenvironment("")
        # add the name of the repair pdb
        if "repairs" in self.params and "pdb" in self.params:
            self.addConfig(
                {"repairpdb": self.params["pdb"] + "_rep" + str(self.params["repairs"])}
            )
        self.addConfig(envs)

    def arg(self, param, override=None):
        if param in self.params:
            return self.params[param]
        elif override != None:
            self.params[param] = override
            return override
        else:
            raise Exception("Missing argument " + param)

    def addConfig(self, config):
        for nm, vl in config.items():
            if nm not in self.params:
                self.params[nm] = vl

    def inputparams(self, argvs, pipeline=False):
        params = {}
        for i in range(1, len(argvs)):
            arg = argvs[i]
            if pipeline:
                params = self.addpipelinetoparams(arg, params)
            else:
                params = self.addlinetoparams(arg, params)
        if "pdb" not in params:
            params["pdb"] = "6vxx"
        if "configfile" not in params:
            params["configfile"] = "../inputs/" + params["pdb"] + "/config.cfg"
        return params

    def spacedparams(self, argvs):
        params = {}
        argvss = argvs[1].split("@")
        for i in range(0, len(argvss)):
            arg = argvss[i]
            params = self.addlinetoparams(arg, params)
        if "pdb" not in params:
            params["pdb"] = "6vxx"
        if "configfile" not in params:
            params["configfile"] = "../inputs/" + params["pdb"] + "/config.cfg"
        return params

    def addlinetoparams(self, arg, params):
        args = []
        if "=" in arg and "@" not in arg:
            args = arg.split("=")
            p, v = args[0], args[1]
            params[p] = v
        return params

    def addpipelinetoparams(self, arg, params):
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

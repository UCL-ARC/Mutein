"""
RSA 5/4/22
------------------------
Class to manage the arguments and precendence of the argument layer
The command line overrides everything
The config file overrides the batch
The batch params ar ethe first level of configutation (script order, dependecies, time, array splits)
"""
import helper as hlp


class ArgumentsX:
    def __init__(self, inputs):
        print("Arguments")
        # combine config and job input params
        iparams = hlp.inputparams(inputs)
        pdb = ""
        if "pdb" in iparams:
            pdb = iparams["pdb"]
        cparams = hlp.configparams(pdb)
        self.params = hlp.mergeparams(cparams, iparams)
        user = self.params["user"]
        user, (foldxe, pythonexe, environment) = hlp.getenvironment(user)
        jobname = self.params["name"]
        input_path, thruput_path, output_path = hlp.get_make_paths(
            pdb, jobname
        )
        self.params["input_path"] = input_path
        self.params["thruput_path"] = thruput_path
        #self.params["interim_path"] = interim_path
        self.params["output_path"] = output_path
        self.params["user"] = user
        self.params["foldxe"] = foldxe
        self.params["pythonexe"] = pythonexe
        self.params["environment"] = environment
        self.params["repairpdb"] = (
            self.params["pdb"] + "_rep" + str(self.params["repairs"])
        )

    def arg(self, param):
        if param in self.params:
            return self.params[param]
        else:
            return ""

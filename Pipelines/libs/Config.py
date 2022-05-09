"""
RSA 5/4/22
------------------------
Class to manage the arguments and precendence of the argument layer
The command line overrides everything
The config file overrides the batch
The batch params ar ethe first level of configutation (script order, dependecies, time, array splits)

This class is an effort to standardise the inputs across pipelines
"""
import yaml


class Config:
    def __init__(self, file_path):
        self.params = {}
        print("Config", file_path)
        with open(file_path, "r") as fr:
            cfgs = yaml.safe_load(fr)
            print("CFG=", cfgs)
            for nm, vl in cfgs.items():
                self.params[nm] = vl

    #############################################################

"""
RSA 23.5.22

Consistent script to run gene level

"""
import os
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)


def prepareGenes(args):
    import Pipelines.geneanalysis.python.mod_datasettogenes as ppla

    ppla.run_pipeline(args)


def preparePdbs(args):
    import Pipelines.geneanalysis.python.mod_datasettopdbs as pplb

    pplb.run_pipeline(args)


def repairPdbs(args):
    print("!!!ERROR irrelevant")


def makeParams(args):
    print("!!!ERROR irrelevant")


def makeVparams(args):
    print("!!!ERROR irrelevant")


def runTasks(args):
    print("!!!ERROR irrelevant")


def runVtasks(args):
    print("!!!ERROR irrelevant")


def aggTasks(args):
    print("!!!ERROR irrelevant")


def aggVtasks(args):
    print("!!!ERROR irrelevant")


def aggGene(args):
    print("!!!ERROR irrelevant")

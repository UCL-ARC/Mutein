"""
RSA 23.5.22

Consistent script to run gene level

"""
import os
import sys
import _helper

def prepareGenes(args):
    import mod_datasettogenes as ppla
    ppla.run_pipeline(args)


def preparePdbs(args):
    import mod_datasettopdbs as pplb

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
    import mod_datasetstitch as ppi
    ppi.run_pipeline(args)

def clean(args):        
    import mod_datasetclean as ppx
    ppx.run_pipeline(args)

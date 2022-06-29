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
    import mod_genetoproteins as ppa
    ppa.run_pipeline(args)


def repairPdbs(args):
    import mod_generepair as ppb
    ppb.run_pipeline(args)


def makeParams(args):
    import mod_geneparams as ppc

    ppc.run_pipeline(args)


def makeVparams(args):
    import mod_genevparams as ppd

    ppd.run_pipeline(args)


def runTasks(args):
    import mod_pdbtask as ppe

    ppe.run_pipeline(args)


def runVtasks(args):
    import mod_pdbvtask as ppf

    ppf.run_pipeline(args)


def aggTasks(args):
    import mod_pdbagg as ppg

    ppg.run_pipeline(args)


def aggVtasks(args):
    import python.mod_pdbvagg as pph

    pph.run_pipeline(args)


def aggGene(args):
    import mod_genestitch as ppi

    ppi.run_pipeline(args)

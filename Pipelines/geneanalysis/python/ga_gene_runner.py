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
    print("!!!ERROR - irrelevant pipeline")

def preparePdbs(args):
    import Pipelines.geneanalysis.python.mod_genetoproteins as ppa
    ppa.run_pipeline(args)
    
def repairPdbs(args):
    import Pipelines.geneanalysis.python.mod_generepair as ppb
    ppb.run_pipeline(args)

def makeParams(args):
    import Pipelines.geneanalysis.python.mod_geneparams as ppc
    ppc.run_pipeline(args)

def makeVparams(args):
    import Pipelines.geneanalysis.python.mod_genevparams as ppd
    ppd.run_pipeline(args)

def runTasks(args):
    import Pipelines.geneanalysis.python.mod_pdbtask as ppe
    ppe.run_pipeline(args)

def runVtasks(args):
    import Pipelines.geneanalysis.python.mod_pdbvtask as ppf
    ppf.run_pipeline(args)

def aggTasks(args):
    import Pipelines.geneanalysis.python.mod_pdbagg as ppg
    ppg.run_pipeline(args)

def aggVtasks(args):
    import Pipelines.geneanalysis.python.mod_pdbvagg as pph
    pph.run_pipeline(args)

def aggGene(args):
    import Pipelines.geneanalysis.python.mod_genestitch as ppi
    ppi.run_pipeline(args)
    
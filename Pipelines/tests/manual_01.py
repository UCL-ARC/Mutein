# test_dummy.py
"""
RSA 11/5/2022
---------------------------
Interim testing of test classes
---------------------------

"""
import os
import sys

import pandas as pd

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)

import SwissModel
import UniProt
import Paths
import PdbRunner

#####################################################
dataset = "play"
gene = "FAT1"
gene = "NOTCH1"


gpth = Paths.Paths("/home/rachel/UCL/github/MuteinData","/home/rachel/UCL/github/Mutein/" + "Pipelines/geneanalysis",dataset=dataset,gene=gene)
# both these searches retun a tuple list of the pdb code and the thruput gene file path, ready for pdb inputs
up = UniProt.UniProt(gene)
df = up.searchForStructures(gpth.gene_outputs,gpth.gene_outpdbs,fragment=50)
accession = up.accession
sm = SwissModel.SwissModel(gene,accession)
dfs = sm.searchForStructures(gpth.gene_outputs,gpth.gene_outpdbs,fragment=50)
dfs.append(df)
vc = pd.concat(dfs, axis=0)
vc.to_csv(gpth.gene_outputs + "Coverage_all.csv",index=False)
print(vc)
pdbs = vc["pdb"].unique()
for pdb in pdbs:
    pth = Paths.Paths("/home/rachel/UCL/github/MuteinData","/home/rachel/UCL/github/Mutein/" + "Pipelines/geneanalysis",dataset=dataset,gene=gene,pdb=pdb)
    prun = PdbRunner.PdbRunner(pdb)
    prun.copyToInput(gpth.gene_outpdbs, pth.pdb_inputs,vc)












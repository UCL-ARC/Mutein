"""
-----------------------------
RSA 15/03/2022
-----------------------------
Adapted from
https://github.com/shorthouse-mrc/COVID_structure/blob/main/Foldx_repair6.py
-----------------------------

This file takes a pdb code (file must be located in the same directory in the format 1xyz.pdb)
It formats the pdb file into a paramater file suitable for foldx PositionScan
The output is in the same directory with the name
scanparams_1xyz.txt
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
from shutil import copyfile
from os.path import exists


# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import Config
import Foldx
import FileDf


def run_pipeline(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print("### FoldX repair job ###")
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset")
    gene = argus.arg("gene")
    pdbcode = argus.arg("pdb").lower()

    pdb_path = Paths.Paths(
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset=dataset,
        gene=gene,
        pdb=pdbcode,
    )

    repair_path = pdb_path.pdb_thruputs + "repair" + str(argus.arg("repairs")) + "/"
    argus.addConfig({"repair_path": repair_path})
    pdb_path.goto_job_dir(repair_path, args, argus.params, "_inputs01")
    ############################################
    pdbfile = pdbcode + ".pdb"
    # Set up files (retain copy of original)
    numRepairs = int(argus.arg("repairs"))
    repair_alreadynames = []
    found = False
    startRep = numRepairs
    base_pdbfile = pdb_path.pdb_inputs + "/" + pdbfile
    while not found and startRep > 0:        
        name = pdb_path.pdb_inputs + "/" + pdbcode + "_rep" + str(startRep) + ".pdb"
        if exists(name):
            found = True
            base_pdbfile = name
        else:
            startRep -= 1

    repairinnames = []
    repairoutnames = []
    for r in range(numRepairs + 1):
        repairinnames.append(pdbcode + "_" + str(r) + ".pdb")
        repairoutnames.append(pdbcode + "_" + str(r) + "_Repair.pdb")

    repairinnames[numRepairs] = pdbcode + "_rep" + str(numRepairs) + ".pdb"
    #### there are 2 files we need in the interim directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
    print(
        "### foldx03: ... copying file",
        base_pdbfile,
        pdb_path.pdb_inputs + "/" + repairinnames[startRep],
        "... ###",
    )

    # first establish which we are currently on

    copyfile(
        base_pdbfile,
        repair_path + repairinnames[startRep],
    )

    # Create Foldx class
    fx_runner = Foldx.Foldx(argus.arg("foldxe"))
    print("Starting repairs from", startRep, "and creating a repair for", numRepairs)
    # Run desired number of repairs
    for r in range(startRep, numRepairs):
        pdb = repairinnames[r]
        output_file = "repair_" + str(r) + ".txt"
        fx_runner.runRepair(pdb, output_file)

        print(
            "### foldx03:  ... copying file",
            argus.arg("repair_path") + repairoutnames[r],
            argus.arg("repair_path") + repairinnames[r + 1],
        )
        copyfile(repairoutnames[r], repairinnames[r + 1])

    # copy the final repaired file to our main interim directory
    print(
        "### ... copying file",
        argus.arg("repair_path") + repairinnames[numRepairs],
        pdb_path.pdb_thruputs + "/" + repairinnames[numRepairs],
    )
    copyfile(
        argus.arg("repair_path") + repairinnames[numRepairs],
        pdb_path.pdb_thruputs + "/" + repairinnames[numRepairs].lower(),
    )
    print(
        "### ... copying file",
        pdb_path.pdb_thruputs + "/" + repairinnames[numRepairs].lower(),
        pdb_path.pdb_inputs + "/" + repairinnames[numRepairs].lower(),
    )
    copyfile(
        pdb_path.pdb_thruputs + "/" + repairinnames[numRepairs].lower(),
        pdb_path.pdb_inputs + "/" + repairinnames[numRepairs].lower(),
    )

    print("### COMPLETED FoldX repair job ###")
    print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)

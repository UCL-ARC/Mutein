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
import sys

import _helper
import Pdb
import Paths
import Arguments
import Config
import Foldx
import FileDf
import AA


def run_pipeline(args):

    ##### INPUTS #############################################
    # The inputs to this function are the pdbfile and the chain id (might optionally consider the positionscan mutation type)
    print("### FoldX repair job ###")
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")    
    data_dir = argus.arg("data_dir")
    dataset = argus.arg("dataset","")
    gene = argus.arg("gene","")
    pdbcode = argus.arg("pdb").lower()

    pdb_path = Paths.Paths(
        data_dir,
        install_dir,
        dataset=dataset,
        gene=gene,
        pdb=pdbcode,
        readonly=False
    )

    repair_path = pdb_path.pdb_thruputs + "repair" + "/"
    repair_from = str(argus.arg("repair_from","x"))
    argus.addConfig({"repair_path": repair_path})
    pdb_path.goto_job_dir(repair_path, args, argus.params, "_inputs01")
    ############################################
    pdbfile = pdbcode + ".pdb"
    # Set up files (retain copy of original)
    numRepairs = int(argus.arg("repairs"))
    repair_alreadynames = []
    found = False
    mostRecentRep = numRepairs
    base_pdbfile = pdb_path.pdb_inputs + "/" + pdbfile
    # find all the existing repair files
    
    while not found:        
        name = repair_path + pdbcode + "_rep" + str(mostRecentRep) + ".pdb"
        name_in = pdb_path.pdb_inputs + pdbcode + "_rep" + str(mostRecentRep) + ".pdb"
        #name = pdb_path.pdb_inputs + "/" + pdbcode + "_rep" + str(mostRecentRep) + ".pdb"
        if exists(name):
            found = True
            base_pdbfile = name                            
            #mostRecentRep += 1        
        elif exists(name_in):
            found = True
            base_pdbfile = name                
            copyfile(name_in,name)        
            #mostRecentRep += 1            
        elif mostRecentRep == 0:
            found = True
        else:
            mostRecentRep -=1
            found = False

    # We are going to start rebuilding as effieciently as possibly from the most recent unless we specifically want to redo
    startRep = mostRecentRep
    if str(repair_from).lower() != "x":        
        startRep = int(repair_from)                    
        base_pdbfile = repair_path + pdbcode + "_" + str(startRep) + ".pdb"

    #repairinnames = []
    #repairoutnames = []
    #repairnextnames = []
    #for r in range(startRep,numRepairs):
    #    repairinnames.append(pdbcode + "_" + str(r) + ".pdb")
    #    repairoutnames.append(pdbcode + "_" + str(r) + "_Repair.pdb")
        #repairnextnames.append(pdbcode + "_" + str(r+1) + ".pdb")

    #repairinnames[numRepairs] = pdbcode + "_rep" + str(numRepairs) + ".pdb"
    #### there are 2 files we need in the interim directory, pdb file rotabase, but rotabase is only needed for foldx4 and NOT needed for foldx5
    if startRep == 0:    
        base_pdbfile = pdb_path.pdb_inputs + "/" + pdbfile
        print(
            "### foldx03: ... copying file",
            base_pdbfile,
            pdb_path.pdb_inputs + "/" + pdbcode + "_rep0.pdb"            
        )
        copyfile(
            base_pdbfile,
            repair_path + pdbcode + "_rep0.pdb"
        )

    # Create Foldx class
    fx_runner = Foldx.Foldx(argus.arg("foldxe"))
    print("Starting repairs from", startRep, "and creating a repair for", numRepairs)
    # Run desired number of repairs
    num_repairs_applied = 0
    lastgoodrepair = numRepairs
    for r in range(startRep, numRepairs):
        pdb = pdbcode + "_rep" + str(r) + ".pdb"
        output_file = "repair_" + str(r) + ".txt"
        return_pdb = pdbcode + "_rep" + str(r) + "_Repair.pdb"
        success = fx_runner.runRepair(pdb, output_file,return_pdb)
        if success:
            
            rename_pdb = pdbcode + "_rep" + str(r+1) + ".pdb"

            print(
                "### foldx03:  ... copying file",
                repair_path + return_pdb,
                repair_path + rename_pdb
            )
            copyfile(repair_path+return_pdb, repair_path+rename_pdb)
            # After every repair check that all the variants exists, if not the last version was our best and we should stop
            filename = pdb_path.pdb_inputs + "params_variants.txt"
            all_variants_exists = True
            pdbobj = Pdb.PdbFile(pdbcode, repair_path + rename_pdb)
            pdbobj.addVariants(filename)
            if pdbobj.existsVariants():            
                all_variants_exists = pdbobj.containsAllVariant()
                
                if not all_variants_exists:
                    lastgoodrepair = r-1
                    break
                num_repairs_applied += 1
            else:
                num_repairs_applied += 1
        
        else:
            lastgoodrepair = r-1
            break

                                                
    # only copy anything if we have been through a loop at all otherwise we have done nothing
    copy_over = True
    if num_repairs_applied == 0:
        # still copy the good file through if it is good
        filename = pdb_path.pdb_inputs + "params_variants.txt"
        last_good_pdb = pdbcode + "_" + str(lastgoodrepair) + ".pdb"
        all_variants_exists = True
        pdbobj = Pdb.PdbFile(pdbcode, last_good_pdb)
        pdbobj.addVariants(filename)
        if pdbobj.existsVariants():            
            all_variants_exists = pdbobj.containsAllVariant()
            if not all_variants_exists:
                copy_over = False    
    if copy_over:
        # copy the final repaired file to our main interim directory
        last_good_pdb = pdbcode + "_rep" + str(lastgoodrepair) + ".pdb"
        save_pdb = pdbcode + "_rep" + str(lastgoodrepair) + ".pdb"
        save_latest_pdb = pdbcode + "_repx.pdb"
        
        try:
            print(
                "### ... copying file",
                repair_path + last_good_pdb,
                pdb_path.pdb_thruputs + save_pdb,
            )
            copyfile(
                repair_path + last_good_pdb,
                pdb_path.pdb_thruputs + "/" + save_pdb
            )
            copyfile(
                repair_path + last_good_pdb,
                pdb_path.pdb_inputs + "/" + save_pdb
            )
            print(
                "### ... copying best file",
                repair_path + last_good_pdb,
                pdb_path.pdb_thruputs + save_latest_pdb,
            )
            copyfile(
                repair_path + last_good_pdb,
                pdb_path.pdb_thruputs + "/" + save_latest_pdb
            )
            copyfile(
                repair_path + last_good_pdb,
                pdb_path.pdb_inputs + "/" + save_latest_pdb
            )
        except:
            print("An error copying the files")
    else:
        print("No repair changes have been made")
    

    print("### COMPLETED FoldX repair job ###")
    #print("MUTEIN SCRIPT ENDED")


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)

"""
RSA 3/5/22
Format calls to foldx
-----------------------------
Calls are either:
repair
posscan
build
-----------------------------
"""

import os
import statistics
import pandas as pd

class Foldx:
    def __init__(self,exe):        
        self.exe = exe        
        # Defaults, thse are explicit here and can be parameterised as and when it is decided to
        self.ionStrength = 0.05
        self.pH = 7
        self.vdwDesign = 2
        self.pdbHydrogens = False

    def makeDefaultsString(self):
        repairB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"
        return repairB
        
    def runRepair(self, pdb, output_file):
        repairA = self.exe + " --command=RepairPDB --pdb="
        repairB = self.makeDefaultsString()
        repairCommand = repairA + pdb + repairB + " > " + output_file
        print("### ... FOLDX:Repair: ", repairCommand)
        os.system(repairCommand)
        print("### ......... FOLDX:Repair: completed")

    def runPosscan(self,pdb,mutation_string):
        foldxcommand = self.exe + " --command=PositionScan "
        foldxcommand += self.makeDefaultsString()                
        foldxcommand += " --pdb=" + pdb
        foldxcommand += " --positions=" + mutation_string
        print("### ... FOLDX:Posscan: ", foldxcommand)
        os.system(foldxcommand)
        print("### ......... FOLDX:Posscan: completed")    

    def runBuild(self,pdb,mutation_string, numRuns)    :
            mut_fl = "individual_list.txt"
            mut_log = "buildmodel.log"
            with open(mut_fl, "w") as fw:
                fw.write(mutation_string + ";")
            # ~/UCL/libs/foldx5/foldx --command=BuildModel --ionStrength=0.05 --pH=7
            #  --water=CRYSTAL --vdwDesign=2 --pdbHydrogens=false --numberOfRuns=15 --mutant-file=mutations.txt
            #  --pdb=6vxx_rep.pdb > buildmodel.log
            foldxcommand = self.exe + " --command=BuildModel "
            foldxcommand += self.makeDefaultsString()                        
            foldxcommand += " --numberOfRuns=" + str(numRuns)
            foldxcommand += " --mutant-file=" + mut_fl
            foldxcommand += " --pdb=" + pdb
            foldxcommand += " > " + mut_log
            print("### ... FOLDX:Build: ", foldxcommand)
            os.system(foldxcommand)
            print("### ......... FOLDX:Build: completed")

    def createBuildCsv(self,pdb,mutation_string,ddg_file,df_file,tag):                    
        ddg_dic = {}
        if os.path.exists(ddg_file):
            with open(ddg_file) as fr:
                jobcontent = fr.readlines()
                energy_list = []
                for linecontents in jobcontent:
                    line = linecontents.split("\t")
                    if line[0].endswith(".pdb"):
                        energy = line[1]
                        energy_list.append(float(energy))
                DDG = statistics.mean(energy_list)
                ddg_dic["mutid"].append(mutation_string.replace(",","_"))
                ddg_dic["ddg"].append(DDG)
                ddg_dic["tag"].append(tag)
                        
            ddg_df = pd.DataFrame.from_dict(ddg_dic)            
            ddg_df.to_csv(df_file, index=False)
            print("### ......... FOLDX:Build_Dataframe: completed",df_file )
        


    

        
        



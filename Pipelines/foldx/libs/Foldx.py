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
from os.path import exists
import statistics
import pandas as pd
import FileDf
import AA


class Foldx:
    def __init__(self, exe):
        self.exe = exe
        # Defaults, thse are explicit here and can be parameterised as and when it is decided to
        self.ionStrength = 0.05
        self.pH = 7
        self.vdwDesign = 2
        self.pdbHydroges = False
        self.AA = AA.AA()
        self.numRuns = 1  # use in BuildModel
        self.outPdb = False

    def makeDefaultsString(self):
        repairB = " --ionStrength=0.05 --pH=7 --vdwDesign=2 --pdbHydrogens=false"
        return repairB

    def runRepair(self, pdb, output_file,return_pdb):
        repairA = self.exe + " --command=RepairPDB --pdb="
        repairB = self.makeDefaultsString()
        repairCommand = repairA + pdb + repairB + " > " + output_file
        print("### ... FOLDX:Repair: ", repairCommand)
        return False
        os.system(repairCommand)        
        if exists(return_pdb):
            print("### ......... FOLDX:Repair: completed")
            return True
        else:
            print("### ......... FOLDX:Repair: failed")
            if exists(output_file):
                with open(output_file, "r") as fr:
                    lines= fr.readlines()
                    for line in lines:
                        print(line.strip())
            return False        

    def runPosscan(self, pdb, mutation_string):
        foldxcommand = self.exe + " --command=PositionScan "
        foldxcommand += self.makeDefaultsString()
        foldxcommand += " --pdb=" + pdb
        foldxcommand += " --positions=" + mutation_string
        print("### ... FOLDX:Posscan: ", foldxcommand)
        os.system(foldxcommand)
        print("### ......... FOLDX:Posscan: completed")

    def runBuild(self, pdb, mutation_strings, tag):
        mut_fl = "individual_list_" + str(tag) + ".txt"
        mut_log = "buildmodel_" + str(tag) + ".log"
        with open(mut_fl, "w") as fw:
            for mutation_string in mutation_strings:
                fw.write(mutation_string + ";\n")
        # ~/UCL/libs/foldx5/foldx --command=BuildModel --ionStrength=0.05 --pH=7
        #  --water=CRYSTAL --vdwDesign=2 --pdbHydrogens=false --numberOfRuns=15 --mutant-file=mutations.txt
        #  --pdb=6vxx_rep.pdb > buildmodel.log
        foldxcommand = self.exe + " --command=BuildModel "
        foldxcommand += self.makeDefaultsString()
        foldxcommand += " --numberOfRuns=" + str(self.numRuns)
        foldxcommand += " --mutant-file=" + mut_fl
        if self.outPdb:
            foldxcommand += " --out-pdb=true"
        else:
            foldxcommand += " --out-pdb=false"
        foldxcommand += " --output-file=" + str(tag)
        foldxcommand += " --pdb=" + pdb
        foldxcommand += " > " + mut_log
        print("### ... FOLDX:Build: ", foldxcommand)
        os.system(foldxcommand)
        print("### ......... FOLDX:Build: completed")

    def createAnaysisCsv(self, pdb, mut_ddg, mut_dic, file_name, score, source):
        posscan_dic = {}
        posscan_dic["source"] = []
        posscan_dic["pdb"] = []
        posscan_dic["score"] = []
        posscan_dic["pdb_mut"] = []
        posscan_dic["pdb_rid"] = []
        posscan_dic["pdb_chain"] = []
        posscan_dic["gene_no"] = []
        posscan_dic["mut_from"] = []
        posscan_dic["mut_to"] = []
        posscan_dic["ddg"] = []
        for pdb_mut, ddg in mut_ddg:
            aaa_from = pdb_mut[:3]
            ch = pdb_mut[3]
            a_to = pdb_mut[-1]
            aaa_to = self.AA.convert(a_to)
            rid = int(pdb_mut[4:-1])
            gene_rid = ""
            if rid in mut_dic:
                gene_rid = mut_dic[rid]
            print(aaa_from, ch, a_to, rid)
            # if the from and to are the same then skip them TODO is this a good decision
            if aaa_from != aaa_to:
                posscan_dic["source"].append(source)
                posscan_dic["pdb"].append(pdb)
                posscan_dic["score"].append(score)
                posscan_dic["pdb_mut"].append(pdb_mut)
                posscan_dic["pdb_rid"].append(rid)
                posscan_dic["pdb_chain"].append(ch)
                posscan_dic["gene_no"].append(gene_rid)
                posscan_dic["mut_from"].append(aaa_from)
                posscan_dic["mut_to"].append(aaa_to)
                posscan_dic["ddg"].append(ddg)

        fdic = FileDf.FileDic(file_name, posscan_dic)
        fdic.saveAsDf()
        print("### ......... FOLDX:Analysis_Dataframe: completed", fdic)

    
    def createBuildCsv(
        # self, pdb, mutation_string, ddg_file, df_file, tag
        self,
        path,
        pdbfile,
        pdb_muts,
        gene_muts,
        coverage,
        ddg_files,
        df_file,
        allAtOnce
    ):
        if allAtOnce:
            self.createBuildCsvAllAtOnce(path,pdbfile,pdb_muts,gene_muts,coverage,ddg_files,df_file)
        else:
            self.createBuildCsvOneByOne(path,pdbfile,pdb_muts,gene_muts,coverage,ddg_files,df_file)
    
    def createBuildCsvAllAtOnce(
        # self, pdb, mutation_string, ddg_file, df_file, tag
        self,
        path,
        pdbfile,
        pdb_muts,
        gene_muts,
        coverage,
        ddg_files,
        df_file        
    ):
        # We have a single results file with nRuns lines per mutation        
        mut_ddg = []
        score, source = 0, "USER"
        if len(coverage["score"]) > 0:
            score = coverage["score"][0]
            source = coverage["source"][0]        
        # there is 1 results file and a list of mutations
        ddg_file = ddg_files[0][0]     
        muts = ddg_files[0][1]
        if exists(ddg_file):  # I want it to error if it is missing
            with open(ddg_file) as fr:
                jobcontent = fr.readlines()
                #skip any lines that are not results
                lines_results = []
                startres= False
                for line in jobcontent:
                    if startres:
                        lines_results.append(line)
                    else:
                        if "Electrostatics" in line:
                            startres = True
                startrow = 0
                for mut_string in muts:                                    
                    energy_list = []
                    for linecontents in lines_results[startrow:startrow+self.numRuns]:
                        line = linecontents.split("\t")
                        if line[0].endswith(".pdb"):
                            energy = line[1]
                            energy_list.append(float(energy))
                    DDG = statistics.mean(energy_list)
                    a_from = mut_string[0]
                    aaa_from = self.AA.convert(a_from)
                    mut_string = aaa_from + mut_string[1:]
                    mut_ddg.append([mut_string, DDG])
                    startrow += self.numRuns
        else:
            raise FileNotFoundError("DDG file does not exist " + ddg_file)

        mut_dic = self.makeMutMapDictionary(pdb_muts, gene_muts, coverage)
        self.createAnaysisCsv(pdbfile, mut_ddg, mut_dic, df_file, score, source)

    
    def createBuildCsvOneByOne(
        # self, pdb, mutation_string, ddg_file, df_file, tag
        self,
        path,
        pdbfile,
        pdb_muts,
        gene_muts,
        coverage,
        ddg_files,
        df_file        
    ):
        # create a list of pdb_muts to ddg
        ### WARNING - there is onyl ONE ddg per build file, unlike posscan where all the ddg are unique
        mut_ddg = []
        score, source = 0, "USER"
        if len(coverage["score"]) > 0:
            score = coverage["score"][0]
            source = coverage["source"][0]        
        for ddg_file, mut_string in ddg_files:
            if exists(ddg_file):  # I want it to error if it is missing
                with open(ddg_file) as fr:
                    jobcontent = fr.readlines()
                    energy_list = []
                    for linecontents in jobcontent:
                        line = linecontents.split("\t")
                        if line[0].endswith(".pdb"):
                            energy = line[1]
                            energy_list.append(float(energy))
                    DDG = statistics.mean(energy_list)
                    a_from = mut_string[0]
                    aaa_from = self.AA.convert(a_from)
                    mut_string = aaa_from + mut_string[1:]
                    mut_ddg.append([mut_string, DDG])
            else:
                raise FileNotFoundError("DDG file does not exist " + ddg_file)

        mut_dic = self.makeMutMapDictionary(pdb_muts, gene_muts, coverage)
        self.createAnaysisCsv(pdbfile, mut_ddg, mut_dic, df_file, score, source)

    def createPosscanCsv(
        self, path, pdbfile, pdb_muts, gene_muts, coverage, ddg_file, outfile_path
    ):
        mut_ddg = []

        fdf = FileDf.FileDf(ddg_file, sep="\t", header=False, cols=["mut", "ddg"])
        df = fdf.openDataFrame()
        score, source = 0, "USER"
        if len(coverage["score"]) > 0:
            score = coverage["score"][0]
            source = coverage["source"][0]
        for i in range(len(df.index)):
            pdb_mut = df["mut"][i]
            ddg = df["ddg"][i]
            mut_ddg.append([pdb_mut, ddg])

        mut_dic = self.makeMutMapDictionary(pdb_muts, gene_muts, coverage)
        self.createAnaysisCsv(pdbfile, mut_ddg, mut_dic, outfile_path, score, source)

    def makeMutMapDictionary(self, pdb_muts, gene_muts, coverage):
        # firs if the pdb_muts' 3rd charavter is a number then convert it to 3 code AA
        for i in range(len(pdb_muts)):
            pdb_mut = pdb_muts[i]
            a = pdb_mut[0]
            ch = pdb_mut[1]
            rid = pdb_mut[2:-1]
            try:
                rid = int(rid)            
                a_from = pdb_mut[0]
                aaa_from = self.AA.convert(a_from)
                pdb_muts[i] = aaa_from + pdb_mut[1:]
            except:
                pass

        mut_dic = {}
        print(pdb_muts)
        if len(pdb_muts) == len(gene_muts):
            for i in range(len(pdb_muts)):  # AA-chain-rid-A
                mut_dic[int(pdb_muts[i][4:-1])] = int(gene_muts[i][2:-1])
                print(int(pdb_muts[i][4:-1]), gene_muts[i][2:-1])
        else:  # then we are going to imply the gene muts from the coverage
            for i in range(len(pdb_muts)):
                pdb_mut = pdb_muts[i]
                ch = pdb_mut[3]
                rid = int(pdb_mut[4:-1])
                print(ch, rid, pdb_mut)
                for i in range(len(coverage.index)):
                    start = int(coverage["pdb_start"][i])
                    end = int(coverage["pdb_end"][i])
                    chain = coverage["chain"][i]
                    gene_start = int(coverage["gene_start"][i])
                    offset = gene_start - start
                    if chain == ch and rid >= start and rid <= end:
                        mut_dic[rid] = int(rid) + int(offset)
                        print(rid, gene_start)

        print(mut_dic)
        return mut_dic

    def getCoveredPdbResidue(self, mut, coverage):
        return ""

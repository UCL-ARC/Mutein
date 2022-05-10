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
import FileDf
import AA


class Foldx:
    def __init__(self, exe):
        self.exe = exe
        # Defaults, thse are explicit here and can be parameterised as and when it is decided to
        self.ionStrength = 0.05
        self.pH = 7
        self.vdwDesign = 2
        self.pdbHydrogens = False
        self.AA = AA.AA()

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

    def runPosscan(self, pdb, mutation_string):
        foldxcommand = self.exe + " --command=PositionScan "
        foldxcommand += self.makeDefaultsString()
        foldxcommand += " --pdb=" + pdb
        foldxcommand += " --positions=" + mutation_string
        print("### ... FOLDX:Posscan: ", foldxcommand)
        os.system(foldxcommand)
        print("### ......... FOLDX:Posscan: completed")

    def runBuild(self, pdb, mutation_string, numRuns):
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

    def createBuildCsv(self, pdb, mutation_string, ddg_file, df_file, tag):
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
                ddg_dic["mutid"].append(mutation_string.replace(",", "_"))
                ddg_dic["ddg"].append(DDG)
                ddg_dic["tag"].append(tag)

            ddg_df = pd.DataFrame.from_dict(ddg_dic)
            ddg_df.to_csv(df_file, index=False)
            print("### ......... FOLDX:Build_Dataframe: completed", df_file)

    def createPosscanCsv(
        self, path, pdbfile, pdb_mut, gene_mut, coverage, outfile_path
    ):
        in_file = path + "PS_" + pdbfile + "_scanning_output.txt"
        print(pdb_mut, gene_mut)
        pdb_muts = pdb_mut.split(",")
        if gene_mut != "x":
            gene_muts = gene_mut.split(",")
        else:
            gene_muts = []
        mut_dic = {}
        print(pdb_muts)
        if len(pdb_muts) == len(gene_muts):
            for i in range(len(pdb_muts)):
                mut_dic[int(pdb_muts[i][2:-1])] = int(gene_muts[i][2:-1])
                print(int(pdb_muts[i][2:-1]), gene_muts[i][2:-1])
        else:  # then we are going to imply the gene muts from the coverage
            for i in range(len(pdb_muts)):
                pdb_mut = pdb_muts[i]
                ch = pdb_mut[1]
                rid = int(pdb_mut[2:-1])
                print(ch, rid, pdb_mut)
                for i in range(len(coverage.index)):
                    start = int(coverage["pdb_start"][i])
                    end = int(coverage["pdb_end"][i])
                    chain = coverage["chain"][i]
                    offset = int(coverage["offset"][i])
                    if chain == ch and rid >= start and rid <= end:
                        mut_dic[rid] = rid + offset
                        print(rid, rid + offset)

        print(mut_dic)

        print(pdb_muts, gene_muts)
        fdf = FileDf.FileDf(in_file, sep="\t", header=False, cols=["mut", "ddg"])
        df = fdf.openDataFrame()
        print(df)
        posscan_dic = {}
        posscan_dic["pdb"] = []
        posscan_dic["pdb_mut"] = []
        posscan_dic["pdb_rid"] = []
        posscan_dic["pdb_chain"] = []
        posscan_dic["gene_no"] = []
        posscan_dic["mut_from"] = []
        posscan_dic["mut_to"] = []
        posscan_dic["ddg"] = []
        for i in range(len(df.index)):
            pdb_mut = df["mut"][i]
            ddg = df["ddg"][i]
            aaa_from = pdb_mut[0:3]
            a_from = self.AA.convert(aaa_from)
            ch = pdb_mut[3]
            a_to = pdb_mut[-1]
            aaa_to = self.AA.convert(a_to)
            rid = int(pdb_mut[4:-1])
            gene_rid = ""
            if rid in mut_dic:
                gene_rid = mut_dic[rid]
            print(aaa_from, ch, a_to, rid)
            posscan_dic["pdb"].append(pdbfile)
            posscan_dic["pdb_mut"].append(pdb_mut)
            posscan_dic["pdb_rid"].append(rid)
            posscan_dic["pdb_chain"].append(ch)
            posscan_dic["gene_no"].append(gene_rid)
            posscan_dic["mut_from"].append(aaa_from)
            posscan_dic["mut_to"].append(aaa_to)
            posscan_dic["ddg"].append(ddg)

        fdic = FileDf.FileDic(outfile_path, posscan_dic)
        fdic.saveAsDf()

    def getCoveredPdbResidue(self, mut, coverage):
        return ""

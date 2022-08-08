"""
RSA 11/5/22
------------------------
Class to manage interactions and data from Swiss Model, the homology structure tool

"""
import os
from urllib.request import urlretrieve
import requests
import genetoprotein
import FileDf
import PdbCoverage
import Pdb

import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)


class SwissModel:
    def __init__(self, gene, accession):
        self.accession = accession
        self.gene = gene

    def searchForStructures(self, gene_path, pdb_path, fragment=-1):
        # An example of the data can be broswed for format in the link below
        # https://swissmodel.expasy.org/repository/uniprot/P46531.json
        url = (
            "https://swissmodel.expasy.org/repository/uniprot/"
            + self.accession.upper()
            + ".json"
        )
        print(url)
        pdb_list = []
        # we will create 2 csv files as a reut, for the homology and for the experimental structures
        fd_hom = FileDf.FileDic(gene_path + "Coverage_sm_hom.csv", {})
        fd_exp = FileDf.FileDic(gene_path + "Coverage_sm_exp.csv", {})

        # sending get request and saving the response as response object
        r = requests.get(url=url)
        # extracting data in json format
        data = r.json()
        res = data["result"]
        struc = res["structures"]
        seq = res["sequence"]

        for s in struc:
            method, reso = "", ""
            frm = str(s["from"])
            to = str(s["to"])
            pdb_file = pdb_path + s["template"] + ".pdb"
            pdb_url = genetoprotein.getPDBLink(s["template"])
            sm_url = s["coordinates"]

            if s["provider"].upper() == "PDB":
                sm_pdb = "smexp_" + s["template"].lower() + "_" + frm + "_" + to
                sm_pdb = sm_pdb.replace(".", "_")
                sm_file = pdb_path + sm_pdb + ".pdb"

                # If it is a pdb structure we walso want the pdb file directly  - WHY? REMOVE THIS
                biopdb = genetoprotein.retrievePdbStructure(
                    pdb_url, s["template"], pdb_file
                )
                smpdb = genetoprotein.retrievePdbStructure(
                    sm_url, s["template"], sm_file
                )
                if biopdb != None:
                    method = biopdb.header["structure_method"]
                    reso = biopdb.header["resolution"]
                pc = PdbCoverage.PdbCoverage(smpdb, seq)
                segments = pc.getCoverage(minfrag=fragment)

                for (
                    chain,
                    residue_num,
                    residue_end,
                    gene_start,
                    gene_end,
                    coverage,
                ) in segments:
                    fd_exp.add("source", "SMEXP")
                    fd_exp.add("gene", self.gene)
                    fd_exp.add("accession", self.accession)
                    fd_exp.add("pdb", sm_pdb)
                    fd_exp.add("method", method)
                    fd_exp.add("resolution", reso)

                    fd_exp.add("chain", chain)
                    fd_exp.add("pdb_start", residue_num)
                    fd_exp.add("pdb_end", residue_end)
                    fd_exp.add("gene_start", gene_start)
                    fd_exp.add("gene_end", gene_end)
                    fd_exp.add("coverage", coverage)
                    fd_exp.add("score", 1)

                if len(segments) > 0:
                    # check it is not CA only
                    onep = Pdb.Pdb(self.gene, sm_pdb, method, reso)
                    onep.segments = segments
                    pdb_fl = Pdb.PdbFile(sm_pdb,sm_file)
                    if not pdb_fl.isCaOnly():
                        pdb_list.append(onep)
                    else:
                        print("CA only pdb file",sm_file)
                    

            elif s["provider"].upper() == "SWISSMODEL":
                sm_pdb = "smhom_" + s["template"].lower() + "_" + frm + "_" + to
                sm_pdb = sm_pdb.replace(".", "_")
                sm_file = pdb_path + sm_pdb + ".pdb"
                smpdb = genetoprotein.retrievePdbStructure(
                    sm_url, s["template"], sm_file
                )
                pc = PdbCoverage.PdbCoverage(smpdb, seq)
                segments = pc.getCoverage(minfrag=fragment)
                similarity = s["similarity"]
                if len(segments) > 0:
                    for (
                        chain,
                        residue_num,
                        residue_end,
                        gene_start,
                        gene_end,
                        coverage,
                    ) in segments:
                        fd_hom.add("source", "SMHOM")
                        fd_hom.add("gene", self.gene)
                        fd_hom.add("accession", self.accession)
                        fd_hom.add("pdb", sm_pdb)
                        fd_hom.add("method", method)
                        fd_hom.add("resolution", reso)

                        fd_hom.add("chain", chain)
                        fd_hom.add("pdb_start", residue_num)
                        fd_hom.add("pdb_end", residue_end)
                        fd_hom.add("gene_start", gene_start)
                        fd_hom.add("gene_end", gene_end)
                        fd_hom.add("coverage", coverage)
                        fd_hom.add("score", similarity)

                    onep = Pdb.Pdb(self.gene, sm_pdb, method, reso)
                    onep.segments = segments
                    # check it is not CA only
                    pdb_fl = Pdb.PdbFile(sm_pdb,sm_file)
                    if not pdb_fl.isCaOnly():
                        pdb_list.append(onep)
                    else:
                        print("CA only pdb file",sm_file)
                    pdb_list.append(onep)

            else:
                print("!!!!", s["template"], s["provider"])

            pdb_url = s["coordinates"]

        df1 = fd_hom.saveAsDf()
        df2 = fd_exp.saveAsDf()
        return [df1, df2], pdb_list

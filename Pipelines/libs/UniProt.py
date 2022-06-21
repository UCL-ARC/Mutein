"""
RSA 11/5/22
------------------------
Class to manage interactions and data from Uniprot, 

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


class UniProt:
    def __init__(self, gene,organism_id):
        self.gene = gene
        self.organism_id = organism_id
        self.accession = genetoprotein.accession_from_bioservices(gene,organism_id)
        seq = genetoprotein.sequence_from_bioservices(self.accession)
        seq_lines = seq.split("\n")
        self.seq = ""
        for s in range(1, len(seq_lines)):
            sl = str(seq_lines[s].strip())
            self.seq += sl

    def searchForStructures(self, gene_path, pdb_path, fragment=-1):
        # using bioservices to get pdbs
        fd_pdb = FileDf.FileDic(gene_path + "Coverage_pdb.csv", {})
        pdb_list = []
        pdbs = genetoprotein.pdbs_from_accession_bioservices(self.accession)
        for pdburl in pdbs:
            pdb = pdburl["pdb"]
            pdb_url = pdburl["path"]
            pdb_file = pdb_path + pdb.lower() + ".pdb"
            biopdb = genetoprotein.retrievePdbStructure(pdb_url, pdb, pdb_file)
            method, reso = "", ""
            if biopdb != None:
                method = biopdb.header["structure_method"]
                reso = biopdb.header["resolution"]
            pc = PdbCoverage.PdbCoverage(biopdb, self.seq)
            segments = pc.getCoverage(minfrag=fragment)
            for (
                chain,
                residue_num,
                residue_end,
                gene_start,
                gene_end,
                coverage,
            ) in segments:
                if pdb[0:2].upper() == "AF":
                    fd_pdb.add("source", "AF")
                else:
                    fd_pdb.add("source", "EXP")
                fd_pdb.add("gene", self.gene)
                fd_pdb.add("accession", self.accession)
                fd_pdb.add("pdb", pdb.lower())
                fd_pdb.add("method", method)
                fd_pdb.add("resolution", reso)

                fd_pdb.add("chain", chain)
                fd_pdb.add("pdb_start", residue_num)
                fd_pdb.add("pdb_end", residue_end)
                fd_pdb.add("gene_start", gene_start)
                fd_pdb.add("gene_end", gene_end)
                fd_pdb.add("coverage", coverage)
                fd_pdb.add(
                    "score", 1
                )  # TODO how to score these, these are real structures not homologous

            if len(segments) > 0:
                onep = Pdb.Pdb(self.gene, pdb.lower(), method, reso)
                onep.segments = segments
                # check it is not CA only
                pdb_fl = Pdb.PdbFile(pdb.lower(),pdb_file)
                if not pdb_fl.isCaOnly():
                    pdb_list.append(onep)
                else:
                    print("CA only pdb file",pdb_file)

        df = fd_pdb.saveAsDf()

        return df, pdb_list

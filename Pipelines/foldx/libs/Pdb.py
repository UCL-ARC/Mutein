"""
RSA 12/4/22
------------------------
Class to manage the data associated with a pdb file specifically for a gene

"""
import os
from os.path import exists
from urllib.request import urlretrieve
import Bio.PDB as bio
import warnings
from Bio import BiopythonWarning

import FileDf
import AA


warnings.simplefilter("ignore", BiopythonWarning)


class Pdb:
    def __init__(
        # self, gene, pdb, chain, segment_start, segment_end, method, resolution,offset
        self,
        gene,
        pdb,
        method,
        resolution,
    ):
        self.gene = gene.upper()
        self.pdbcode = pdb
        self.segments = []
        # self.chain = chain.upper()
        # self.segment_starts = []  # [int(segment_start)]
        # self.segment_ends = []  # [int(segment_end)]
        # self.segment_offsets = []  # [int(offset)]
        # self.segment_chains = []  # [chain]
        # self.method = method.lower()
        # self.resolution = resolution

    """
    def addSegment(self, seg_start, seg_end, seg_off, seg_chain):
        start_in = seg_start not in self.segment_starts
        end_in = seg_end not in self.segment_ends
        chain_in = seg_chain not in self.segment_chains
        if start_in and end_in:
            self.segment_starts.append(int(seg_start))
            self.segment_ends.append(int(seg_end))
            self.segment_offsets.append(int(seg_off))
            self.segment_chains.append(seg_chain)

    """

    def matchesResidue(self, residue):
        # for s in range(len(self.segment_starts)):
        for (
            chain,
            residue_num,
            residue_end,
            gene_start,
            gene_end,
            coverage,
        ) in self.segments:
            if int(residue) >= residue_num and int(residue) <= residue_end:
                return True
        return False

    def getSegments(self):
        segs = ""
        for (
            chain,
            residue_num,
            residue_end,
            gene_start,
            gene_end,
            coverage,
        ) in self.segments:
            segs += (
                chain
                + ":"
                + residue_num
                + ":"
                + residue_end
                + ":"
                + gene_start
                + ":"
                + gene_end
                + ":"
                + coverage
            )
            # segs += seg_start + ":" + seg_end + ":" + seg_off + ":" + seg_chain + " "
        return segs

    def getMethod(self):
        if "nmr" in self.method:
            return "nmr"
        elif "x-ray" in self.method:
            return "x-ray"
        elif "electron microscopy" in self.method:
            return "e-m"
        elif "unknown" in self.method and "AF" in self.pdbcode:
            return "alphafold"
        else:
            return self.method

    ### These functions retrieve pdb structures from the web ###

    def downloadPdb(self, file_path, location="pdb"):  # pdb, pdbe, af
        if self.getMethod() == "alphafold":
            web_path = self.getAlphaFoldLink()
        elif location == "pdbe":
            web_path = self.getPDBeLink()
        else:
            web_path = self.getPDBLink()

        return self.retrievePdbStructure(
            web_path, file_path + "/" + self.pdbcode.lower() + ".pdb"
        )

    def retrievePdbStructure(self, url, path_name):
        if self.retrieveFile(url, path_name):
            parser = bio.PDBParser()
            struc = parser.get_structure(self.pdbcode, path_name)
            return struc
        return None

    def removePdbStructure(self, url, pdb, path_name):
        if os.path.exists(path_name):
            os.remove(path_name)

    def retrieveFile(self, url, path_name):
        if not os.path.exists(path_name):
            try:
                print(url)
                urlretrieve(url, path_name)
                return True
            except:
                print("...!!! No data for", url)
                return False
        return True

    ### These functions define the web links where the structures can be downloaded from ###
    def getPDBeLink(self):
        return (
            "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb"
            + self.pdbcode.lower()
            + ".ent"
        )

    def getPDBLink(self):
        return "https://files.rcsb.org/download/" + self.pdbcode.upper() + ".pdb"

    def getAlphaFoldLink(self):
        af_path = "https://alphafold.ebi.ac.uk/files/" + self.pdbcode + ".pdb"
        return af_path

class PdbFile:
    def __init__(self,pdbcode,filepath):                    
        self.pdbcode = pdbcode
        self.filepath = filepath
        self.lines = []

        if exists(self.filepath):            
            with open(self.filepath) as rf:
                self.lines = rf.readlines()                
        
        self.variants = []

    def addVariants(self,variant_file):
        if exists(variant_file):
            fdf = FileDf.FileDf(variant_file, sep=" ")
            csv = fdf.openDataFrame()                        
            pdb_muts = csv["pdb_mut"]
            #all_muts = []
            for pmut in pdb_muts:
                muts = pmut.split(",")
                for mut in muts:
                    self.variants.append(mut)

    def existsVariants(self):
        return len(self.variants) > 0
    
    def isCaOnly(self):
        if self.pdbcode == "1pet":
            bk = "yes"
        ca = True
        for line in self.lines:
            line = line.strip()
            if len(line) > 15:
                atm = line[0:4]
                if atm.upper() == "ATOM":
                    if "CA" not in line:
                        ca = False                
        return ca

    def containsAllVariant(self):
        aa_conv = AA.AA()
        contains_all = True
        for mut in self.variants:
            a = mut[0]
            aaa = aa_conv.convert(a)
            rid = mut[2:-1]
            chain = mut[1]
            string_match = aaa + " " + chain
            if int(rid) < 10:
                string_match += "   "
            elif int(rid) < 100:
                string_match += "  "
            elif int(rid) < 1000:
                string_match += " "
            string_match += rid
            this_variant = False
            for line in self.lines:
                line = line.strip()
                if string_match in line:
                    this_variant = True
            if not this_variant:
                contains_all = False
        return contains_all


"""
RSA 12/4/22
------------------------
Class to manage the data associated with a pdb file specifically for a gene

"""
import os
from urllib.request import urlretrieve
import Bio.PDB as bio
import warnings

warnings.filterwarnings("ignore")  # sadly because of the annoying bioython warnings


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
        # self.chain = chain.upper()
        self.segment_starts = []  # [int(segment_start)]
        self.segment_ends = []  # [int(segment_end)]
        self.segment_offsets = []  # [int(offset)]
        self.segment_chains = []  # [chain]
        self.method = method.lower()
        self.resolution = resolution

    def addSegment(self, seg_start, seg_end, seg_off, seg_chain):
        start_in = seg_start not in self.segment_starts
        end_in = seg_end not in self.segment_ends
        chain_in = seg_chain not in self.segment_chains
        if start_in and end_in:
            self.segment_starts.append(int(seg_start))
            self.segment_ends.append(int(seg_end))
            self.segment_offsets.append(int(seg_off))
            self.segment_chains.append(seg_chain)

    def matchesResidue(self, residue):
        for s in range(len(self.segment_starts)):
            seg_start = int(self.segment_starts[s])
            seg_end = int(self.segment_ends[s])
            seg_off = int(self.segment_offsets[s])
            seg_chain = self.segment_chains[s]
            if int(residue) >= seg_start and int(residue) <= seg_end:
                return True
        return False

    def getSegments(self):
        segs = ""
        for s in range(len(self.segment_starts)):
            seg_start = str(self.segment_starts[s])
            seg_end = str(self.segment_ends[s])
            seg_off = str(self.segment_offsets[s])
            seg_chain = str(self.segment_chains[s])
            segs += seg_start + ":" + seg_end + ":" + seg_off + ":" + seg_chain + " "
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

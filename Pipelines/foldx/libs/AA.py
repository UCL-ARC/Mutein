"""
RSA 4/5/22
This just performs conversion between types of amino acid code

"""


class AA:
    def __init__(self):
        self.three = {
            "ALA": "A",
            "CYS": "C",
            "ASP": "D",
            "GLU": "E",
            "PHE": "F",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LYS": "K",
            "LEU": "L",
            "MET": "M",
            "ASN": "N",
            "PRO": "P",
            "GLN": "Q",
            "ARG": "R",
            "SER": "S",
            "THR": "T",
            "VAL": "V",
            "TRP": "W",
            "TYR": "Y",
            "H1S": "o",
            "H2S": "e",
        }
        self.one = {
            "A": "ALA",
            "C": "CYS",
            "D": "ASP",
            "E": "GLU",
            "F": "PHE",
            "G": "GLY",
            "H": "HIS",
            "I": "ILE",
            "K": "LYS",
            "L": "LEU",
            "M": "MET",
            "N": "ASN",
            "P": "PRO",
            "Q": "GLN",
            "R": "ARG",
            "S": "SER",
            "T": "THR",
            "V": "VAL",
            "W": "TRP",
            "Y": "TYR",
            "o": "H1S",
            "e": "H2S",
        }

    def convert(self, aa):
        if len(aa) == 1:
            if aa in self.one:
                return self.one[aa]
        elif len(aa) == 3:
            if aa in self.three:
                return self.three[aa]
        return aa

"""
RSA 30/6/22

This only needs to be done once, we go through the 38 transcript and save the chromosomes, genes and the start and end codon positions

"""
import os
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(lib_path)
import RefGenome


import pandas as pd
rg = RefGenome.RefGenome()
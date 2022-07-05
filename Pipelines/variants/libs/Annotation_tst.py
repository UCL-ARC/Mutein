
"""
RSA 30/6/22

class to test RefGenome
"""

import pandas as pd
import Annotation

ann_path = "/home/rachel/UCL/github/MuteinData/data_genome/gencode.v40.annotation.gff3"

anno = Annotation.Annotation(ann_path)

genes = []
genes.append("ANG")
genes.append("APOE")
genes.append("APP")
genes.append("NOTCH1")

anno.getCdsRegions(genes)

        
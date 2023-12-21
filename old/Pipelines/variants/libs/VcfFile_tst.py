
"""
RSA 30/6/22

class to test RefGenome
"""

import pandas as pd
import VcfFile

hardcoded_vcf_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/keogh2018.vcf"
vf = VcfFile.VcfFile(hardcoded_vcf_path)
df = vf.getVariants()
print(df)




        
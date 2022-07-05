
"""
RSA 30/6/22

class to test RefGenome
"""

import pandas as pd
import Fasta

fasta_path = "/home/rachel/UCL/github/MuteinData/data_genome/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
fst = Fasta.Fasta(fasta_path)
# FROM shearwater
print(fst.getSeq(1,7961913,7961913)=="G")
print(fst.getSeq(1,7965215,7965215)=="C")
print(fst.getSeq(1,7969434,7969434)=="T")
print(fst.getSeq(1,7969449,7969449)=="G")
print(fst.getSeq(1,7971065,7971065)=="G")
print(fst.getSeq(10,43102747,43102747)=="A")
print(fst.getSeq(10,43105241,43105241)=="A")
print(fst.getSeq(10,43106301,43106301)=="T")
print(fst.getSeq(11,32417455,32417455)=="T")
print(fst.getSeq(11,32417492,32417492)=="T")
print(fst.getSeq(11,32427871,32427871)=="C")
print(fst.getSeq(11,32427874,32427874)=="C")
print(fst.getSeq(6,152061259,152061259)=="C")
print(fst.getSeq(6,152061285,152061285)=="A")
print(fst.getSeq(6,152098960,152098960)=="G")
print(fst.getSeq(7,92615108,92615108)=="C")
print(fst.getSeq(7,92618019,92618019)=="A")
print(fst.getSeq(7,92622935,92622935)=="A")
print(fst.getSeq(7,95301079,95301079)=="A")
print(fst.getSeq(7,95308134,95308134)=="T")
print(fst.getSeq(7,95411704,95411704)=="G")
print(fst.getSeq(7,95412230,95412230)=="G")
print(fst.getSeq(9,136506452,136506452)=="C")
print(fst.getSeq(9,136507052,136507052)=="G")
print(fst.getSeq(9,136507194,136507194)=="T")
print(fst.getSeq(9,136508211,136508211)=="T")



        
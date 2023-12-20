
"""
RSA 30/5/22
Test the results
"""
import pandas as pd
import Analysis

# First create the dataframes
#ddg_background = "data/ddg_background.csv"
ddg_background = "data/7w7g_background.csv"
ddg_var_bm = "data/ddg_var_bm.csv"
ddg_var_ps = "data/ddg_var_ps.csv"

mode = "PDB_BACK"

df_back = pd.read_csv(ddg_background)
df_varbm = pd.read_csv(ddg_var_bm)
df_varps = pd.read_csv(ddg_var_ps)

ana = Analysis.Analysis(df_back)
if mode=="GENE":
    ana.createPdbSummary()
    ana.createPdbSummary()
    ana.show()    

elif mode=="PDB_BACK":
    ana.createDdgBackRid()
    ana.createDdgBackFromTo()
    ana.histAllBackground()
    ana.show()


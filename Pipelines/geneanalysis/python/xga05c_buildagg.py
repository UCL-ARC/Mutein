"""
-----------------------------
RSA 17/03/22
-----------------------------

This aggregates the outputs from split positionscans into 1 file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import statistics
import pandas as pd
from shutil import copyfile

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Paths
import Arguments
import Config


def run_pipeline07(args):
    ##### INPUTS #############################################
    print("### Foldx variant aggregate ###")
    ##############################################
    argus = Arguments.Arguments(args)
    install_dir = argus.arg("install_dir")
    sys.path.append(install_dir)
    sys.path.append(install_dir + "/Pipelines")
    sys.path.append(install_dir + "/Pipelines/libs")
    data_dir = argus.arg("data_dir")
    pdbcode = argus.arg("pdb").lower()
    pdb_path = Paths.Paths(
        "pdb",
        data_dir,
        install_dir + "Pipelines/geneanalysis",
        dataset="",
        gene="",
        pdb=pdbcode,
    )
    pdb_config = Config.Config(pdb_path.pdb_inputs + "/config.yml")
    argus.addConfig(pdb_config.params)
    work_path = pdb_path.pdb_thruputs + "vagg/"
    argus.addConfig({"work_path": work_path})
    pdb_path.goto_job_dir(argus.arg("work_path"), args, argus.params, "_inputs07")

    params_file = pdb_path.pdb_thruputs + "variant_params.txt"
    variant_dirs = []
    with open(params_file) as fr:
        paramscontent = fr.readlines()
        for pc in paramscontent:
            pcs = pc.strip().split(" ")
            vnam = pcs[2].replace(",", "_")
            vdir = pcs[3]
            variant_dirs.append([vdir, vnam])

    ddg_dic = {"mutid": [], "ddg": [], "tag": []}
    for vd, vn in variant_dirs:
        print(vd, vn)
        ddg_file = (
            pdb_path.pdb_thruputs + vd + "/Dif_" + argus.arg("repairpdb") + ".fxout"
        )  # the pdb repaired files are always pdbcode_rep
        print(ddg_file)
        if os.path.exists(ddg_file):
            with open(ddg_file) as fr:
                jobcontent = fr.readlines()
                energy_list = []
                for linecontents in jobcontent:
                    line = linecontents.split("\t")
                    if line[0].endswith(".pdb"):
                        energy = line[1]
                        energy_list.append(float(energy))
                DDG = statistics.mean(energy_list)
                ddg_dic["mutid"].append(vn)
                ddg_dic["ddg"].append(DDG)
                ddg_dic["tag"].append(vd)

    # Make a dataframe
    import pandas as pd

    ddg_df = pd.DataFrame.from_dict(ddg_dic)
    df_file = argus.arg("repairpdb") + "_variants_ddg_dataframe.csv"
    ddg_df.to_csv(pdb_path.pdb_outputs + df_file, index=False)
    print("saved dataframe to", pdb_path.pdb_outputs + df_file)

    # And save something visual as a starting point for some analysis
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 9))
    fig.suptitle(
        argus.arg("pdb") + " variant mutations\nddg <-1=stabilising >2.5=destabilising"
    )
    sns.set_color_codes("pastel")
    # first plt
    sns.barplot(x="ddg", y="tag", data=ddg_df, color="g", ax=ax1)
    ax1.set_ylabel("")
    # second plt
    sns.histplot(data=ddg_df, x="ddg", palette="tab20", ax=ax2, bins=20)
    ax2.set_ylabel("")
    sns.despine(left=True, bottom=True)
    plt.rcParams["axes.labelsize"] = 25
    plot_file = argus.arg("repairpdb") + "_variant_plot.png"
    plt.savefig(pdb_path.pdb_outputs + plot_file)
    print("saved plot to", pdb_path.pdb_outputs + plot_file)

    print("### COMPLETED FoldX aggregate job ###")
    print("MUTEIN SCRIPT ENDED")


#####################################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline07"](sys.argv)

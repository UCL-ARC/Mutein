"""
RSA 11/4/22
------------------------
Class to manage the path hierarchy of the pipelines

This class is an effort to standardise the inputs across pipelines

If any restructuring is wanted, it should be possible to do it all here and feed through automatically.
"""
# import helper as hlp
import os


class PathsXXX:
    def __init__(self, data_dir, pipe_path, dataset="", gene="", pdb=""):
        # Depending on the levels we will ensure the correct paths are present for inputs and ouputs
        self.data_dir = data_dir
        if self.data_dir[-1] != "/":
            self.data_dir += "/"
        print("Initialising PATH CLASS with datadir=", data_dir)
        self.pipeline_path = pipe_path

        # there are levels of ids as many to 1 is slowly whittled down by dataset, gene, pdb and perhaps method (eg ddg or other)
        # [dataset] 1--* [gene] 1--* [pdb] 1--* [method]

        # Create the dataset paths
        level = "pdb"
        if pdb == "":
            level = "gene"
            if gene == "":
                level = "dataset"
        dataset = dataset.replace(".","_")
        gene = gene.replace(".","_")
        pdb = pdb.replace(".","_")

        dataset = dataset.lower()
        gene = gene.lower()
        pdb = pdb.lower()

        if "dataset" == level:
            # clean up first?
            self.get_make_paths_vcf(dataset)
        # Create the geneprot paths
        elif "gene" == level:
            # clean up first?
            self.get_make_paths_geneprot(dataset, gene)
        elif "pdb" == level:
            # clean up first?
            self.get_make_paths_pdb(dataset, gene, pdb)

        # Create the foldx paths (which may be unique for a method)

    def get_make_paths_vcf(self, dataset):
        paths_dic = {}
        # The path structure is relative to the script file (it will be easy to change to a given one)

        in_root_path = "inputs"
        print(in_root_path)
        self.add_remove_path_levels(0, dir=in_root_path)
        in_root_path = "inputs/datasets"
        self.add_remove_path_levels(0, dir=in_root_path)

        out_root_path = "outputs"
        self.add_remove_path_levels(0, dir=out_root_path)
        out_root_path = "outputs/datasets"
        self.add_remove_path_levels(0, dir=out_root_path)

        thru_root_path = "thruputs"
        self.add_remove_path_levels(0, dir=thru_root_path)
        thru_root_path = "thruputs/datasets"
        self.add_remove_path_levels(0, dir=thru_root_path)

        self.dataset_inputs = self.add_remove_path_levels(
            0, dir=in_root_path + "/" + dataset + "/"
        )
        self.dataset_outputs = self.add_remove_path_levels(
            0, dir=out_root_path + "/" + dataset + "/"
        )

    def get_make_paths_geneprot(self, dataset, gene):
        paths_dic = {}
        # The path structure is relative to the script file (it will be easy to change to a given one)
        dir_script_path = os.path.dirname(os.path.realpath(__file__))
        dir_path = self.add_remove_path_levels(2)

        in_root_path = "inputs"
        self.add_remove_path_levels(0, dir=in_root_path)
        in_root_path = "inputs/genes"
        self.add_remove_path_levels(0, dir=in_root_path)

        out_root_path = "outputs"
        self.add_remove_path_levels(0, dir=out_root_path)
        out_root_path = "outputs/genes"
        self.add_remove_path_levels(0, dir=out_root_path)

        thru_root_path = "thruputs"
        self.add_remove_path_levels(0, dir=thru_root_path)
        thru_root_path = "thruputs/genes"
        self.add_remove_path_levels(0, dir=thru_root_path)

        append = ""
        if dataset != "":
            append = (dataset + "_").lower()

        self.gene_inputs = self.add_remove_path_levels(
            0,
            dir=in_root_path + "/" + append + gene + "/",
            must_exist=False,
        )
        self.gene_outputs = self.add_remove_path_levels(
            0, dir=out_root_path + "/" + append + gene + "/"
        )
        self.add_remove_path_levels(0, dir=thru_root_path + "/" + append + gene)
        self.gene_outpdbs = self.add_remove_path_levels(
            0, dir=thru_root_path + "/" + append + gene + "/pdbs/"
        )
        # inputs and outputs directories are overkill here

    def get_make_paths_pdb(self, dataset, gene, pdbcode):
        # The path structure is relative to the script file (it will be easy to change to a given one)
        in_root_path = "inputs"
        self.add_remove_path_levels(0, dir=in_root_path)
        in_root_path = "inputs/pdbs"
        self.add_remove_path_levels(0, dir=in_root_path)

        out_root_path = "outputs"
        self.add_remove_path_levels(0, dir=out_root_path)
        out_root_path = "outputs/pdbs"
        self.add_remove_path_levels(0, dir=out_root_path)

        thru_root_path = "thruputs"
        self.add_remove_path_levels(0, dir=thru_root_path)
        thru_root_path = "thruputs/pdbs"
        self.add_remove_path_levels(0, dir=thru_root_path)

        append = ""
        if dataset != "":
            append = (dataset + "_").lower()
        if gene != "":
            append += gene + "_"

        self.pdb_inputs = self.add_remove_path_levels(
            0, dir=in_root_path + "/" + append + pdbcode.lower() + "/"
        )
        self.pdb_thruputs = self.add_remove_path_levels(
            0, dir=thru_root_path + "/" + append + pdbcode.lower() + "/"
        )
        self.pdb_outputs = self.add_remove_path_levels(
            0, dir=out_root_path + "/" + append + pdbcode.lower() + "/"
        )

    def add_remove_path_levels(self, levels=0, dir="", must_exist=False):
        dirs = self.data_dir.split("/")
        if levels > 0:
            dirs = dirs[: -1 * levels]
        retpath = "/".join(dirs) + dir
        if os.path.exists(retpath):
            return retpath
        elif must_exist:
            raise Exception("The folder is missing: " + retpath)
        else:
            try:
                os.mkdir(retpath)
            except:
                pass #this could be a rare case of 2 nodes creating it at the same time
        return retpath

    def goto_job_dir(self, dir_path, args, params, name):
        if not os.path.exists(dir_path):
            try:
                os.mkdir(dir_path)
            except:
                pass #this could be a rare case of 2 nodes creating it at the same time
        os.chdir(dir_path)
        inputs_file = name + ".log"
        with open(inputs_file, "w") as fw:
            for arg in args:
                fw.write(str(arg) + " ")
            fw.write("\n")
            for cfg, val in params.items():
                fw.write(str(cfg) + "=" + str(val) + "\n")

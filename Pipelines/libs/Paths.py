"""
RSA 11/4/22
------------------------
Class to manage the path hierarchy of the pipelines

This class is an effort to standardise the inputs across pipelines

If any restructuring is wanted, it should be possible to do it all here and feed through automatically.
"""
# import helper as hlp
import os


class Paths:
    def __init__(self, data_dir, pipe_path, dataset="", gene="", pdb=""):
        # Depending on the levels we will ensure the correct paths are present for inputs and ouputs
        self.dataset_inputs = ""
        self.dataset_thruputs = ""
        self.dataset_outputs = ""

        self.gene_inputs = ""
        self.gene_thruputs = ""
        self.gene_outputs = ""

        self.pdb_inputs = ""
        self.pdb_thruputs = ""
        self.pdb_outputs = ""
                
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
            self.get_make_paths_dataset(dataset)
        # Create the geneprot paths
        elif "gene" == level:
            # clean up first?
            self.get_make_paths_gene(dataset, gene)
        elif "pdb" == level:
            # clean up first?
            self.get_make_paths_pdb(dataset, gene, pdb)

        # Create the foldx paths (which may be unique for a method)

    def get_make_paths_dataset(self, dataset):        
        # The path structure is relative to the script file (it will be easy to change to a given one)
        basepath = "dataset_" + dataset + "/"
        self.dataset_inputs = "dataset_" + dataset + ""
        print(self.dataset_inputs)
        self.dataset_inputs = self.add_remove_path_levels(0, dir=basepath)
                
        self.dataset_outputs = basepath + "results"
        self.dataset_outputs = self.add_remove_path_levels(0, dir=self.dataset_outputs)
        
        self.dataset_thruputs = basepath + "thruputs"
        self.dataset_thruputs = self.add_remove_path_levels(0, dir=self.dataset_thruputs)
                

    def get_make_paths_gene(self, dataset, gene):        
        # The path structure is relative to the script file (it will be easy to change to a given one)
        basepath = ""
        if dataset != "":
            self.get_make_paths_dataset(dataset)
            basepath = "dataset_" + dataset
        
        if basepath == "" or basepath[-1] != "/":
            basepath += "/"
                
        self.gene_inputs = basepath + "gene_" + gene + ""
        self.gene_inputs = self.add_remove_path_levels(0, dir=self.gene_inputs)
        basepath = basepath + "gene_" + gene + "/"
        
        self.gene_outputs = basepath + "results"
        self.gene_outputs = self.add_remove_path_levels(0, dir=self.gene_outputs)

        self.gene_thruputs = basepath + "thruputs"
        self.gene_thruputs = self.add_remove_path_levels(0, dir=self.gene_thruputs)

        self.gene_outpdbs = basepath + "thruputs/pdbs"
        self.gene_outpdbs = self.add_remove_path_levels(0, dir=self.gene_outpdbs)
                        
    def get_make_paths_pdb(self, dataset, gene, pdbcode):
        # The path structure is relative to the script file (it will be easy to change to a given one)            
        basepath = ""
        if dataset != "":
            self.get_make_paths_dataset(dataset)
            basepath = basepath = "dataset_" + dataset + "/"
        if gene != "":
            self.get_make_paths_gene(dataset,gene)
            basepath += "gene_" + gene + "/"
                
        self.pdb_inputs = basepath + "pdb_" + pdbcode + ""
        self.pdb_inputs = self.add_remove_path_levels(0, dir=self.pdb_inputs)
        basepath = basepath + "pdb_" + pdbcode + "/"
        
        self.pdb_outputs = basepath + "results"
        self.pdb_outputs = self.add_remove_path_levels(0, dir=self.pdb_outputs)

        self.pdb_thruputs = basepath + "thruputs"
        self.pdb_thruputs = self.add_remove_path_levels(0, dir=self.pdb_thruputs)

    def add_remove_path_levels(self, levels=0, dir="", must_exist=False):
        dirs = self.data_dir.split("/")
        if levels > 0:
            dirs = dirs[: -1 * levels]
        retpath = "/".join(dirs) + dir
        if retpath[-1] == "/":
            retpath[:-1]
        if os.path.exists(retpath):
            return retpath + "/"
        elif must_exist:
            raise Exception("The folder is missing: " + retpath)
        else:
            try:                
                os.mkdir(retpath)            
            except:
                pass #this could be a rare case of 2 nodes creating it at the same time
        return retpath + "/"

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

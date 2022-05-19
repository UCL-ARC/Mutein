"""
RSA 3/5/22
Automatically create batches of this format
-----------------------------
module load python3/recommended
python pipeline_qsubber.py batch_tst01.yml qsub notch NOTCH1 ""

"""
import os
import subprocess
import SubRunner


class BatchMaker:
    def __init__(self, script_file, yaml_file):
        self.batches = []
        self.batches.append(
            "### SCRIPT AUTOMATICALLY GENERATED BY MUTEIN PIPELINE (UCL-ARC 2022) ###"
        )
        self.batches.append("run=qsub")
        self.batches.append("install_dir=$1")
        self.batches.append("script=${install_dir}Pipelines/" + script_file)
        self.batches.append("config=${install_dir}Pipelines/" + yaml_file)
        self.batches.append('echo "EXE PATH=$install_dir"')
        self.batches.append('echo "CURRENT=$PWD"')
        #self.batches.append("module load python3/recommended")       
        self.batches.append(
            'echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"'
        )

    def addBatch(self, dataset, gene):
        # python ${script} $install_dir $PWD ${config} $run notch NOTCH1 1pb5
        line = (
            "python ${script} $install_dir $PWD ${config} $run "
            + dataset
            + " "
            + gene            
        )
        self.batches.append(line)

    def printBatchScript(self, file_path, sym_link=""):
        #self.batches.append('echo "MUTEIN SCRIPT ENDED"')
        print("BATCH created in", file_path)
        with open(file_path, "w") as fw:
            for line in self.batches:
                fw.write(line + "\n")
        self.changeMod(file_path)
        # if sym_link !="":
        #    self.createSymLink(file_path,sym_link)
        #    self.changeMod(sym_link)
        return file_path
    
    def printBatches(self,file_path,batches, script_dir):
        print("BATCH created in", file_path)
        with open(file_path, "w") as fw:
            for line in batches:
                fw.write(line + " " + script_dir + "\n")
        self.changeMod(file_path)

    def createSymLink(self, file_path, sym_path):
        # ln -s ~/code/notes/notes ~/bin/notes
        args = []
        args.append("ln")
        args.append("-s")
        args.append(file_path)
        args.append(sym_path)
        print("### subprocess():", args)
        process = subprocess.Popen(
            args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

    def changeMod(self, file_path):
        # ln -s ~/code/notes/notes ~/bin/notes
        args = []
        args.append("chmod")
        args.append("+x")
        args.append(file_path)
        print("### subprocess():", args)
        process = subprocess.Popen(
            args=args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

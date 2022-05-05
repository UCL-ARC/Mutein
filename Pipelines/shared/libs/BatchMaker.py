"""
RSA 3/5/22
Automatically create batches of this format
-----------------------------
module load python3/recommended
python pipeline_qsubber.py batch_tst01.yml qsub notch NOTCH1 ""

"""



class BatchMaker:
    def __init__(self):        
        self.batches = []
        self.batches.append("module load python3/recommended")

    def addBatch(self,script_file,yaml_file,dataset,gene,pdb,chain):
        line = "python " + script_file + " " + yaml_file + " qsub " + dataset + " " + gene + " " + pdb + " " + chain
        self.batches.append(line)

    def printBatchScript(self,file_path):
        with open(file_path, "w") as fw:        
            for line in self.batches:
                fw.write(line + "\n")
            

    

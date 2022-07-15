"""
RSA: 19/5/22
This CI test script runs from the pdb level
It enables debugging of the scripts as if run from a batch

"""
import sys, os
# import from the shared library in Mutein/Pipelines/shared/lib
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
python_path = "/".join(dirs) + "/python"
lib_path = "/".join(dirs) + "/libs"
sys.path.append(python_path)
sys.path.append(lib_path)
import Paths

def addpath(inputs):        
    print(sys.path)
    inputs += "@install_dir=/home/rachel/UCL/github/Mutein/"
    inputs += "@data_dir=/home/rachel/UCL/github/MuteinData/"
    return inputs

######################################################################
### INPUTS
#dataset=""
#gene=""
#pdb="1pb5"

dataset="mouse"
gene="nup214"
pdb=""#smhom_4ll7_3_b_805_882"

repairs="4"
repair_from = "x"
split=100
vsplit=10
task=3
runs = "e"
missing="N"
variant="*"
""" 
    if "a" in runs:        print("Mutein: Preparing genes")        
    if "b" in runs:        print("Mutein: Preparing pdbs")        
    if "c" in runs:        print("Mutein: Repairing pdbs")            
    if "d" in runs:        print("Mutein: Making background param file")        
    if "e" in runs:        print("Mutein: Making variant param file")        
    if "f" in runs:        print("Mutein: Running background tasks")            
    if "g" in runs:        print("Mutein: Running variant tasks")            
    if "h" in runs:        print("Mutein: Aggregating background tasks")            
    if "i" in runs:        print("Mutein: Aggregating variant tasks")            
    if "j" in runs:        print("Mutein: Aggregating gene tasks")            
    if "k" in runs:        print("Mutein: Aggregating dataset")        
    
    if "x" in runs:        print("Mutein: Cleaning data aggressively and packaging batch")        
"""

# which steps of the pipeline to run
def runPPL(inputs):
    inputs = addpath(inputs)
    import ga__runner as ppl
    args = ["", inputs]
    ppl.run_pipeline(args)

runPPL("runs="+runs+"@dataset="+dataset+"@gene="+gene+"@pdb="+pdb+"@repairs="+str(repairs)+"@repair_from="+str(repair_from)+"@task=" +str(task)+"@chunk="+str(split)+"@variant=*@vchunk="+str(vsplit)+"@missing="+missing+"@variant="+variant)


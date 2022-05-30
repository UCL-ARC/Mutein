
import tkinter as tk
from functools import partial

global username
global password
global server
global pipeline_dir
global working_dir
global install_dir
global data_dir

username = "ucbxxxx"
password = "xyz"
server = "myriad"
pipeline_dir = "/home/ucbtlcr/Mutein/Pipelines/geneanalysis/"
working_dir = "/home/ucbtlcr/Scratch/workspace/"
install_dir = "/home/ucbtlcr/Mutein/"
data_dir = "/home/ucbtlcr/MuteinData/"


def getSshShell():
    global username
    global password
    global server
    import spur
    try:
        shell = spur.SshShell(hostname= server + ".rc.ucl.ac.uk", username=username, password=password)
        return shell
    except:
        return None

def RunClean():    
        shell = getSshShell()   
        global pipeline_dir
        global install_dir      
        scriptname = pipeline_dir + "REM_errors.sh"
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "CLEAN",install_dir])
        return result.output

def RunRerun():    
        shell = getSshShell()   
        global pipeline_dir
        global install_dir      
        scriptname = pipeline_dir + "REM_errors.sh"
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "RERUNERROR",install_dir])
        return result.output

def RunCleanAll():    
        shell = getSshShell()                
        global pipeline_dir
        global install_dir      
        scriptname = pipeline_dir + "REM_errors.sh"
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "ALL",install_dir])
        return result.output

def RunView(jobid):        
        shell = getSshShell()
        global pipeline_dir
        global install_dir      
        scriptname = pipeline_dir + "REM_view.sh"
        result = shell.run(["chmod","+x",scriptname])        
        result = shell.run([scriptname, jobid,install_dir])
        return result.output

def SubmitJob(command,dataset,gene,pdb):        
    global pipeline_dir
    global working_dir
    global install_dir
    global data_dir
    if dataset == "":
        dataset = "x"
    if gene == "":
        gene = "x"
    if pdb == "":
        pdb = "x"
    #mode, pattern, WorkDir, DataDir, InstallDir, PipelineDir
    shell = getSshShell()
    scriptname = ""
    if command == "DATASET_PDBS":
        scriptname = pipeline_dir + "REMOTE_foldx_dataset_pdbs.sh"
    elif command == "GENE_PDBS":
        scriptname = pipeline_dir + "REMOTE_foldx_gene_pdbs.sh"
    elif command == "GENE_REP":
        scriptname = pipeline_dir + "REMOTE_foldx_gene_rep.sh"
    elif command == "GENE_TASKS":
        scriptname = pipeline_dir + "REMOTE_foldx_gene_tasks.sh"
    elif command == "GENE_SPLIT":
        scriptname = pipeline_dir + "REMOTE_foldx_gene_split.sh"
    elif command == "GENE_AGG":
        scriptname = pipeline_dir + "REMOTE_foldx_gene_agg.sh"
    elif command == "PDB_REP":
        scriptname = pipeline_dir + "REMOTE_foldx_pdb_rep.sh"
    elif command == "PDB_SPLIT":
        scriptname = pipeline_dir + "REMOTE_foldx_pdb_split.sh"
    elif command == "PDB_TASKS":
        scriptname = pipeline_dir + "REMOTE_foldx_pdb_tasks.sh"
    elif command == "PDB_UNTASKS":
        scriptname = pipeline_dir + "REMOTE_foldx_pdb_untasks.sh"
    elif command == "PDB_AGG":
        scriptname = pipeline_dir + "REMOTE_foldx_pdb_agg.sh"        
    if scriptname != "":
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, dataset,gene,pdb,working_dir,data_dir,install_dir,pipeline_dir])
        return result.output
        #return scriptname
    else:
        return "Script not yet implemented " + command

def RunScript(mode,pattern):    
    global pipeline_dir
    global working_dir
    global install_dir
    global data_dir
    #mode, pattern, WorkDir, DataDir, InstallDir, PipelineDir
    scriptname = pipeline_dir + "REMOTE.sh"
    shell = getSshShell()                
    result = shell.run(["chmod","+x",scriptname])
    result = shell.run([scriptname, mode,pattern,working_dir,data_dir,install_dir,pipeline_dir])
    return result.output
    

def make_paths(txtBox,txtUser,txtPwd,txtServer):
    global username
    global password
    global server
    global pipeline_dir
    global working_dir
    global install_dir
    global data_dir
    username = txtUser.get().strip().lower()
    password = txtPwd.get().strip()
    server = txtServer.get().strip()
    pipeline_dir = "/home/"+username+"/Mutein/Pipelines/geneanalysis/"
    working_dir = "/home/"+username+"/Scratch/workspace/"
    install_dir = "/home/"+username+"/Mutein/"
    data_dir = "/home/"+username+"/MuteinData/"
    
    path_show = "USER =\t\t\t"+username
    ssh_path = server + ".rc.ucl.ac.uk"
    path_show += "\nSSH =\t\t\t" + ssh_path    
    path_show += "\nWorkDir =\t\t\t" + working_dir 
    path_show += "\nDataDir =\t\t\t" + data_dir    
    path_show += "\nInstallDir =\t\t\t" + install_dir    
    path_show += "\nPipelineDir =\t\t\t" + pipeline_dir           
        
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', path_show)

def verify_connection(txtBox):
    shell = getSshShell()
    if shell == "None":
        result = "!!!FAILED!!!\nCheck your inputs, connection and VPN"
    else:
        try:
            result = shell.run(["ls"]).output.decode(encoding="utf-8")    
            if len(result) > 1:
                result = "SUCCESS\nThe home contents are\n\n" + str(result)
            else:
                result = "!!!FAILED!!!\nCheck your inputs, connection and VPN"
        except:
            result = "!!!FAILED!!!\nCheck your inputs, connection and VPN"
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', result)
    
            
class tabConfig:

    def __init__(self,parent):
        self.parent = parent
        self.txtBox = None
        self.txtJob = None
                                                
    def createTab(self, clr,clra):                
        ### layout ###
        
        lhs = tk.Label(self.parent,bg=clra,text="")
        lhs.pack(padx=5, pady=15, side=tk.LEFT)
        rhs = tk.Label(self.parent,bg=clra,text="")
        rhs.pack(padx=5, pady=15, side=tk.LEFT)
        
        header = tk.Label(rhs,bg=clra,text="User config for Mutein")
        header.grid(row=0,column=0,columnspan=2,sticky = 'EW')
        
        self.txtBox = tk.Text(rhs,height=25,width=80,borderwidth=5)        
        self.txtBox.grid(row=1,column=0,rowspan=10)        
        self.txtBox.insert('end', "")
        scrollb = tk.Scrollbar(rhs,command=self.txtBox.yview)
        scrollb.grid(row=1, column=1, sticky='nsew',rowspan=10)
        self.txtBox['yscrollcommand'] = scrollb.set
                        
        lbl_user = tk.Label(lhs,bg=clra,text="User name",width=20)
        lbl_pwd = tk.Label(lhs,bg=clra,text="Password")
        lbl_server = tk.Label(lhs,bg=clra,text="Server",width=20)        
        lbl_empty = tk.Label(lhs,bg=clra,text="",width=20)        
        text_user = tk.Entry(lhs,width=20,borderwidth=5, relief="sunken")        
        #text_pwd = tk.Text(lhs,height=1,width=20,borderwidth=5, relief="sunken")                
        text_pwd = tk.Entry(lhs,show="*",width=20,borderwidth=5, relief="sunken")                        
        text_server = tk.Entry(lhs,width=20,borderwidth=5, relief="sunken")

        text_user.insert(0, username)
        text_pwd.insert(0, password)
        text_server.insert(0, server)

        btnG = tk.Button(lhs, text="Regenerate Paths", command=partial(make_paths,self.txtBox,text_user,text_pwd,text_server),borderwidth=5, relief="groove",width=20)
        btnP = tk.Button(lhs, text="Verify Connection", command=partial(verify_connection,self.txtBox),borderwidth=5, relief="groove",width=20)        
                
        btnG.grid(row=6,column=1)        
        btnP.grid(row=7,column=1)        
                                
                
        lbl_user.grid(row=1,column=0)                
        lbl_pwd.grid(row=2,column=0)        
        lbl_server.grid(row=3,column=0)        
        lbl_empty.grid(row=4,column=0)                
        text_user.grid(row=1,column=1)                
        text_pwd.grid(row=2,column=1)        
        text_server.grid(row=3,column=1)        
        




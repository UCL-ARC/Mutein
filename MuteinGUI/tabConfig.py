
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
# Directories on the server
pipeline_trunc = "/Scratch/Mutein/Pipelines/foldx/"
scripts_trunc = "/Scratch/Mutein/Pipelines/foldx/scripts/"
working_trunc = "/Scratch/workspace/"
install_trunc = "/Scratch/Mutein/"
data_trunc = "/Scratch/MuteinData/Foldx/"
#data_trunc = "/MuteinData/"
# global dir variables
pipeline_dir = "/home/"+username+pipeline_trunc
scripts_dir = "/home/"+username+scripts_trunc
working_dir = "/home/"+username+working_trunc
install_dir = "/home/"+username+install_trunc
data_dir = "/home/"+username+data_trunc
# scripts on the server
remove_script = "REM_errors.sh"
view_script = "REM_view.sh"
remote_script = "REMOTE.sh"
scriptname_ds_pdb = "REMOTE_foldx_dataset_pdbs.sh"
scriptname_ds_agg = "REMOTE_foldx_dataset_agg.sh"    
scriptname_gene_pdb = "REMOTE_foldx_gene_pdbs.sh"
scriptname_gene_rep = "REMOTE_foldx_gene_rep.sh"    
scriptname_gene_unrep = "REMOTE_foldx_gene_unrep.sh"    
scriptname_gene_tasks = "REMOTE_foldx_gene_tasks.sh"    
scriptname_gene_untasks = "REMOTE_foldx_gene_untasks.sh" 
scriptname_gene_split = "REMOTE_foldx_gene_split.sh"    
scriptname_gene_agg = "REMOTE_foldx_gene_agg.sh"    
scriptname_pdb_pdb = "REMOTE_foldx_pdb_pdb.sh"    
scriptname_pdb_rep = "REMOTE_foldx_pdb_rep.sh"    
scriptname_pdb_split = "REMOTE_foldx_pdb_split.sh"    
scriptname_pdb_tasks = "REMOTE_foldx_pdb_tasks.sh"    
scriptname_pdb_un = "REMOTE_foldx_pdb_untasks.sh"    
scriptname_pdb_agg = "REMOTE_foldx_pdb_agg.sh"        
scriptname_ds_clean = "REMOTE_foldx_dataset_clean.sh"        
scriptname_gene_clean = "REMOTE_foldx_gene_clean.sh"        
scriptname_pdb_clean = "REMOTE_foldx_pdb_clean.sh"        


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

def sendSshShell(script,cmds):
    global username
    global password
    global server
    import spur
    try:
        shell = spur.SshShell(hostname= server + ".rc.ucl.ac.uk", username=username, password=password)
        result = shell.run(["chmod","+x",script])
        cms = [script]
        for cmd in cmds:
            cms.append(cmd)        
        print("SHELL",cms)
        result = shell.run(cms)        
        return result
    except:
        return None

def RunClean():    
        #shell = getSshShell()   
        global pipeline_dir
        global scripts_dir
        global install_dir      
        scriptname = scripts_dir + remove_script
        #result = shell.run(["chmod","+x",scriptname])
        #result = shell.run([scriptname, "CLEAN",install_dir])
        #return result.output
        return sendSshShell(scriptname,["CLEAN",install_dir]).output

def RunRerun():    
        shell = getSshShell()   
        global pipeline_dir
        global scripts_dir
        global install_dir      
        scriptname = scripts_dir + remove_script
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "RERUNERROR",install_dir])
        return result.output

def RunRerunO():    
        shell = getSshShell()   
        global pipeline_dir
        global scripts_dir
        global install_dir      
        scriptname = scripts_dir + remove_script
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "RERUNORPHANS",install_dir])
        return result.output

def RunCleanAll():    
        shell = getSshShell()                
        global pipeline_dir
        global scripts_dir
        global install_dir      
        scriptname = scripts_dir + remove_script
        result = shell.run(["chmod","+x",scriptname])
        result = shell.run([scriptname, "ALL",install_dir])
        return result.output

def RunView(jobid):        
        shell = getSshShell()
        global pipeline_dir
        global scripts_dir
        global install_dir      
        scriptname = scripts_dir + view_script
        result = shell.run(["chmod","+x",scriptname])        
        result = shell.run([scriptname, jobid,install_dir])
        return result.output

def SubmitJob(command,dataset,gene,pdb):        
    global pipeline_dir
    global scripts_dir
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
    genes = gene.split(",")
    ret = ""
    for gene in genes:
        shell = getSshShell()
        scriptname = ""
        if command == "DATASET_PDBS":
            scriptname = scripts_dir + scriptname_ds_pdb
        elif command == "DS_CLEAN":
            scriptname = scripts_dir + scriptname_ds_clean
        elif command == "GENE_PDBS":
            scriptname = scripts_dir + scriptname_gene_pdb
        elif command == "GENE_REP":
            scriptname = scripts_dir + scriptname_gene_rep
        elif command == "GENE_UNREP":
            scriptname = scripts_dir + scriptname_gene_unrep
        elif command == "GENE_TASKS":
            scriptname = scripts_dir + scriptname_gene_tasks
        elif command == "GENE_UNTASKS":
            scriptname = scripts_dir + scriptname_gene_untasks
        elif command == "GENE_SPLIT":
            scriptname = scripts_dir + scriptname_gene_split
        elif command == "GENE_AGG":
            scriptname = scripts_dir + scriptname_gene_agg
        elif command == "DS_AGG":
            scriptname = scripts_dir + scriptname_ds_agg
        elif command == "GENE_CLEAN":
            scriptname = scripts_dir + scriptname_gene_clean
        elif command == "PDB_PDB":
            scriptname = scripts_dir + scriptname_pdb_pdb
        elif command == "PDB_REP":
            scriptname = scripts_dir + scriptname_pdb_rep
        elif command == "PDB_SPLIT":
            scriptname = scripts_dir + scriptname_pdb_split
        elif command == "PDB_TASKS":
            scriptname = scripts_dir + scriptname_pdb_tasks
        elif command == "PDB_UNTASKS":
            scriptname = scripts_dir + scriptname_pdb_un
        elif command == "PDB_AGG":
            scriptname = scripts_dir + scriptname_pdb_agg
        elif command == "PDB_CLEAN":
            scriptname = scripts_dir + scriptname_pdb_clean
        if scriptname != "":
            result = shell.run(["chmod","+x",scriptname])
            result = shell.run([scriptname, dataset,gene,pdb,working_dir,data_dir,install_dir,pipeline_dir])
            ret += (result.output).decode('utf-8')
            #return scriptname
        else:
            ret += "Script not yet implemented " + command
    return ret

def RunScript(mode,pattern):    
    global pipeline_dir
    global working_dir
    global install_dir
    global data_dir
    #mode, pattern, WorkDir, DataDir, InstallDir, PipelineDir
    scriptname = scripts_dir + remote_script    
    result = sendSshShell(scriptname,[mode,pattern,working_dir,data_dir,install_dir,pipeline_dir])
    return result.output
    
def make_paths(txtBox,txtUser,txtPwd,txtServer):
    global username
    global password
    global server
    global pipeline_dir
    global scripts_dir
    global working_dir
    global install_dir
    global data_dir
    username = txtUser.get().strip().lower()
    password = txtPwd.get().strip()
    server = txtServer.get().strip()
    pipeline_dir = "/home/"+username+pipeline_trunc
    scripts_dir = "/home/"+username+scripts_trunc
    working_dir = "/home/"+username+working_trunc
    install_dir = "/home/"+username+install_trunc
    data_dir = "/home/"+username+data_trunc
    
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
        




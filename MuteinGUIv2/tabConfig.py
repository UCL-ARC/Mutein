
import tkinter as tk
from functools import partial
import SshComms as sc

                        
class tabConfig:

    def __init__(self,app,parent,txtView):
        self.app = app
        self.parent = parent
        self.txtBox = txtView
        self.txtJob = None

        self.username = "ucbxxxx"
        self.password = "xyz"
        self.server = "myriad"
        self.pipeline_dir = "/home/ucbtlcr/Mutein/Pipelines/geneanalysis/"
        self.working_dir = "/home/ucbtlcr/Scratch/workspace/"
        self.install_dir = "/home/ucbtlcr/Mutein/"
        self.data_dir = "/home/ucbtlcr/MuteinData/"
                                                    
    def createTab(self, clr,clra):                
        ### layout ###
        
        lhs = tk.Label(self.parent,bg=clra,text="")
        lhs.pack(padx=5, pady=15, side=tk.LEFT)
                                                                
        lbl_user = tk.Label(lhs,bg=clra,text="User name",width=20)
        lbl_pwd = tk.Label(lhs,bg=clra,text="Password")
        lbl_server = tk.Label(lhs,bg=clra,text="Server",width=20)        
        lbl_empty = tk.Label(lhs,bg=clra,text="",width=20)        
        text_user = tk.Entry(lhs,width=20,borderwidth=5, relief="sunken")        
        #text_pwd = tk.Text(lhs,height=1,width=20,borderwidth=5, relief="sunken")                
        text_pwd = tk.Entry(lhs,show="*",width=20,borderwidth=5, relief="sunken")                        
        text_server = tk.Entry(lhs,width=20,borderwidth=5, relief="sunken")

        text_user.insert(0, self.username)
        text_pwd.insert(0, self.password)
        text_server.insert(0, self.server)

        btnG = tk.Button(lhs, text="Regenerate Paths", command=partial(self.make_paths,text_user,text_pwd,text_server),borderwidth=5, relief="groove",width=20)
        btnP = tk.Button(lhs, text="Verify Connection", command=partial(self.verify_connection),borderwidth=5, relief="groove",width=20)        
                
        btnG.grid(row=6,column=1)        
        btnP.grid(row=7,column=1)        
                                
                
        lbl_user.grid(row=1,column=0)                
        lbl_pwd.grid(row=2,column=0)        
        lbl_server.grid(row=3,column=0)        
        lbl_empty.grid(row=4,column=0)                
        text_user.grid(row=1,column=1)                
        text_pwd.grid(row=2,column=1)        
        text_server.grid(row=3,column=1)     

    def runSshShell(self,command):        
        import SshComms as sc         
        download_thread = sc.AsyncSsh()
        download_thread.setDetails(self.username,self.password,self.server,command)
        download_thread.start()
        self.app.monitor(download_thread)        

    def RunClean(self):    
            shell = getSshShell()   
            global pipeline_dir
            global install_dir      
            scriptname = pipeline_dir + "REM_errors.sh"
            result = shell.run(["chmod","+x",scriptname])
            result = shell.run([scriptname, "CLEAN",install_dir])
            return result.output

    def RunRerun(self):    
            shell = getSshShell()   
            global pipeline_dir
            global install_dir      
            scriptname = pipeline_dir + "REM_errors.sh"
            result = shell.run(["chmod","+x",scriptname])
            result = shell.run([scriptname, "RERUNERROR",install_dir])
            return result.output

    def RunCleanAll(self):    
            shell = getSshShell()                
            global pipeline_dir
            global install_dir      
            scriptname = pipeline_dir + "REM_errors.sh"
            result = shell.run(["chmod","+x",scriptname])
            result = shell.run([scriptname, "ALL",install_dir])
            return result.output

    def RunView(self,jobid):        
            shell = getSshShell()
            global pipeline_dir
            global install_dir      
            scriptname = pipeline_dir + "REM_view.sh"
            result = shell.run(["chmod","+x",scriptname])        
            result = shell.run([scriptname, jobid,install_dir])
            return result.output

    def SubmitJob(self,command,dataset,gene,pdb):        
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

    def RunScript(self,mode,pattern):    
        global pipeline_dir
        global working_dir
        global install_dir
        global data_dir
        #mode, pattern, WorkDir, DataDir, InstallDir, PipelineDir
        scriptname = pipeline_dir + "REMOTE.sh"
        runSshShell(["chmod","+x",scriptname])                
        runSshShell([scriptname, mode,pattern,working_dir,data_dir,install_dir,pipeline_dir])                
            
    def make_paths(self,txtUser,txtPwd,txtServer):
        
        self.username = txtUser.get().strip().lower()
        self.password = txtPwd.get().strip()
        self.server = txtServer.get().strip()

        self.pipeline_dir = "/home/"+self.username+"/Mutein/Pipelines/geneanalysis/"
        self.working_dir = "/home/"+self.username+"/Scratch/workspace/"
        self.install_dir = "/home/"+self.username+"/Mutein/"
        self.data_dir = "/home/"+self.username+"/MuteinData/"
        
        path_show = "USER =\t\t\t"+self.username
        self.ssh_path = self.server + ".rc.ucl.ac.uk"
        path_show += "\nSSH =\t\t\t" + self.ssh_path    
        path_show += "\nWorkDir =\t\t\t" + self.working_dir 
        path_show += "\nDataDir =\t\t\t" + self.data_dir    
        path_show += "\nInstallDir =\t\t\t" + self.install_dir    
        path_show += "\nPipelineDir =\t\t\t" + self.pipeline_dir           
            
        self.txtBox.delete(1.0,tk.END)
        self.txtBox.insert('end', path_show)

    def verify_connection(self):
        self.runSshShell(["ls"])   
            






import tkinter as tk
from functools import partial
import tabConfig as rs
           
def show_genes(txtBox,txtDataset,txtGene,txtPdb):    
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("GENES",dataset+":"+gene+":"+pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def show_pdbs(txtBox,txtDataset,txtGene,txtPdb):    
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("PDBS",dataset+":"+gene+":"+pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def show_status(txtBox,txtDataset,txtGene,txtPdb):    
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("PDB",dataset+":"+gene+":"+pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_pdbs(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = ""
    pdb = ""
    ret = rs.SubmitJob("DATASET_PDBS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_rep(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_REP",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_unrep(txtBox,txtDataset,txtGene,txtPdb):        
    run_ds_missing_pdbs(txtBox,txtDataset,txtGene,txtPdb)
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_UNREP",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_split(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_SPLIT",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_tasks(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_TASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_untasks(txtBox,txtDataset,txtGene,txtPdb):        
    run_ds_renew_tasks(txtBox,txtDataset,txtGene,txtPdb)
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_UNTASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_agg(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = "ALL"
    pdb = ""
    ret = rs.SubmitJob("GENE_AGG",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_clean(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = ""
    pdb = ""
    ret = rs.SubmitJob("DS_CLEAN",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_pdbs(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_PDBS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_rep(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_REP",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_unrep(txtBox,txtDataset,txtGene,txtPdb):        
    run_gene_missing_pdbs(txtBox,txtDataset,txtGene,txtPdb)
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_UNREP",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_split(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_SPLIT",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_tasks(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_TASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_untasks(txtBox,txtDataset,txtGene,txtPdb):        
    run_gene_renew_tasks(txtBox,txtDataset,txtGene,txtPdb)
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("GENE_UNTASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_agg(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_AGG",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_clean(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = ""
    ret = rs.SubmitJob("GENE_CLEAN",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_pdb(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()    
    ret = rs.SubmitJob("PDB_PDB",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_rep(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_REP",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_split(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_SPLIT",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_missing_pdbs(txtBox,txtDataset,txtGene,txtPdb):            
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()    
    ret = rs.RunScript("GENEINCOMPLETEPDB",dataset+":"+gene+":x")
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_missing_pdbs(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = "ALL"
    ret = rs.RunScript("DSINCOMPLETEPDB",dataset+":x:x")
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_renew_tasks(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("PDBINCOMPLETE",dataset+":"+gene+":"+pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_gene_renew_tasks(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("GENEINCOMPLETE",dataset+":"+gene+":x")
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_ds_renew_tasks(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.RunScript("DSINCOMPLETE",dataset+":x:x")
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_tasks(txtBox,txtDataset,txtGene,txtPdb):            
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_TASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_untasks(txtBox,txtDataset,txtGene,txtPdb):        
    run_pdb_renew_tasks(txtBox,txtDataset,txtGene,txtPdb)
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_UNTASKS",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_agg(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_AGG",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)
def run_pdb_clean(txtBox,txtDataset,txtGene,txtPdb):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    ret = rs.SubmitJob("PDB_CLEAN",dataset,gene,pdb)
    txtBox.delete(1.0,tk.END)
    txtBox.insert('end', ret)

        
class tabPipeline:

    def __init__(self,parent):
        self.parent = parent
        self.txtBox = None
        self.txtJob = None
                                                
    def createTab(self,clr,clra):                
        ### layout ###
        
        lhs = tk.Label(self.parent,bg=clra,text="")
        lhs.pack(padx=5, pady=15, side=tk.LEFT)
        rhs = tk.Label(self.parent,bg=clra,text="")
        rhs.pack(padx=5, pady=15, side=tk.LEFT)
        
        header = tk.Label(lhs,bg=clra,text="Pipeline Control for Mutein")
        header.grid(row=0,column=0,columnspan=2,sticky = 'EW')
        
        self.txtBox = tk.Text(lhs,height=40,width=100,borderwidth=5)        
        self.txtBox.grid(row=1,column=0,rowspan=10)        
        self.txtBox.insert('end', "")
        scrollb = tk.Scrollbar(lhs,command=self.txtBox.yview)
        scrollb.grid(row=1, column=1, sticky='nsew',rowspan=10)
        self.txtBox['yscrollcommand'] = scrollb.set

                
        header2 = tk.Label(rhs,bg="silver",text="Run batch")
        lblDataset = tk.Label(rhs,bg=clra,text="dataset",width=20)
        lblGene = tk.Label(rhs,bg=clra,text="gene",width=20)
        lblPdb = tk.Label(rhs,bg=clra,text="pdb",width=20)

        lblWarning = tk.Label(rhs,bg=clra,text="Warning, only submit new batches to the pipeline if you are really sure!",width=60)
        lblMissing = tk.Label(rhs,bg=clra,text="This submits just missing tasks to the pipeline",width=60)
        lblClean = tk.Label(rhs,bg=clra,text="Take care, this deletes ALL INTERIM FILES",width=60)

        text_dataset = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")        
        text_gene = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")        
        text_pdb = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")

        text_dataset.insert(0, "sheartop")
        text_gene.insert(0, "NOTCH1")
        text_pdb.insert(0, "1pb5")
                
        clrAlert = "red"
        clrAlertAlert = "tomato"
        clrGoForIt = "darkseagreen"
        clrMaybe = "lightsalmon"

        btnG = tk.Button(rhs, text="View dataset status", command=partial(show_genes,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg=clrGoForIt)
        btnP = tk.Button(rhs, text="View gene status", command=partial(show_pdbs,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg=clrGoForIt)
        btnS = tk.Button(rhs, text="View pdb status", command=partial(show_status,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg=clrGoForIt)
                
        # Pipeline Dataset
        btnSubmit_dspdbs = tk.Button(rhs, text="Submit pdb prepare", command=partial(run_ds_pdbs,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg=clrMaybe)
        btnSubmit_dsrep = tk.Button(rhs, text="Submit repair", command=partial(run_ds_rep,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_dsmrep = tk.Button(rhs, text="Submit missing repair", command=partial(run_ds_unrep,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_dssplit = tk.Button(rhs, text="Submit splits prepare", command=partial(run_ds_split,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_dstasks = tk.Button(rhs, text="Submit tasks", command=partial(run_ds_tasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)
        btnSubmit_duntasks = tk.Button(rhs, text="Submit missing tasks", command=partial(run_ds_untasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)
        btnSubmit_dsagg = tk.Button(rhs, text="Submit aggregation", command=partial(run_ds_agg,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_dsclean = tk.Button(rhs, text="Submit clean", command=partial(run_ds_clean,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg="orange")
        
        # Pipeline Gene
        btnSubmit_gppdbs = tk.Button(rhs, text="Submit pdb prepare", command=partial(run_gene_pdbs,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg=clrMaybe)
        btnSubmit_grep = tk.Button(rhs, text="Submit repair", command=partial(run_gene_rep,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_gmrep = tk.Button(rhs, text="Submit missing repair", command=partial(run_gene_unrep,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_gprep = tk.Button(rhs, text="Submit splits prepare", command=partial(run_gene_split,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_gtasks = tk.Button(rhs, text="Submit tasks", command=partial(run_gene_tasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)
        btnSubmit_guntasks = tk.Button(rhs, text="Submit missing tasks", command=partial(run_gene_untasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)
        btnSubmit_gagg = tk.Button(rhs, text="Submit aggregation", command=partial(run_gene_agg,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_gclean = tk.Button(rhs, text="Submit clean", command=partial(run_gene_clean,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg="orange")

        # Pipeline PDB
        btnSubmit_pdbpdb = tk.Button(rhs, text="Download pdb", command=partial(run_pdb_pdb,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_rep = tk.Button(rhs, text="Submit repair", command=partial(run_pdb_rep,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlertAlert)
        btnSubmit_prep = tk.Button(rhs, text="Submit splits prepare", command=partial(run_pdb_split,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_ptasks = tk.Button(rhs, text="Submit all tasks", command=partial(run_pdb_tasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)        
        btnSubmit_puntasks = tk.Button(rhs, text="Submit missing tasks", command=partial(run_pdb_untasks,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrAlert)        
        btnSubmit_pagg = tk.Button(rhs, text="Submit aggregation", command=partial(run_pdb_agg,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg=clrMaybe)
        btnSubmit_pclean = tk.Button(rhs, text="Submit clean", command=partial(run_pdb_clean,self.txtBox,text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20, bg="orange")
        
        header2.grid(row=0,column=0,columnspan=6,sticky = 'EW')        
        lblDataset.grid(row=1,column=0, padx=2, pady=2)    
        lblGene.grid(row=1,column=1, padx=2, pady=2)
        lblPdb.grid(row=1,column=2, padx=2, pady=2)

        btnG.grid(row=2,column=0, padx=2, pady=2)
        btnP.grid(row=2,column=1, padx=2, pady=2)
        btnS.grid(row=2,column=2, padx=2, pady=2)

        text_dataset.grid(row=3,column=0, padx=2, pady=2)      
        text_gene.grid(row=3,column=1, padx=10, pady=2)
        text_pdb.grid(row=3,column=2, padx=2, pady=2)

        lblWarning.grid(row=4,column=0, padx=2, pady=2,columnspan=3)
                        
        btnSubmit_dspdbs.grid(row=5,column=0, padx=2, pady=2)
        btnSubmit_dsrep.grid(row=6,column=0, padx=2, pady=2)        
        btnSubmit_dssplit.grid(row=7,column=0, padx=2, pady=2)
        btnSubmit_dstasks.grid(row=8,column=0, padx=2, pady=2)                
        btnSubmit_dsagg.grid(row=9,column=0, padx=2, pady=2)
         
        btnSubmit_gppdbs.grid(row=5,column=1, padx=2, pady=2)
        btnSubmit_grep.grid(row=6,column=1, padx=2, pady=2)
        btnSubmit_gprep.grid(row=7,column=1, padx=2, pady=2)
        btnSubmit_gtasks.grid(row=8,column=1, padx=2, pady=2)
        btnSubmit_gagg.grid(row=9,column=1, padx=2, pady=2)
                
        btnSubmit_pdbpdb.grid(row=5,column=2, padx=2, pady=2)
        btnSubmit_rep.grid(row=6,column=2, padx=2, pady=2)
        btnSubmit_prep.grid(row=7,column=2, padx=2, pady=2)
        btnSubmit_ptasks.grid(row=8,column=2, padx=2, pady=2)                
        btnSubmit_pagg.grid(row=9,column=2, padx=2, pady=2)

        lblMissing.grid(row=10,column=0, padx=2, pady=2,columnspan=3)
        btnSubmit_dsmrep.grid(row=11,column=0, padx=2, pady=2)
        btnSubmit_duntasks.grid(row=12,column=0, padx=2, pady=2)
        btnSubmit_gmrep.grid(row=11,column=1, padx=2, pady=2)
        btnSubmit_guntasks.grid(row=12,column=1, padx=2, pady=2)
        btnSubmit_puntasks.grid(row=12,column=2, padx=2, pady=2)
        # the cleaning which deletes all interim results
        lblClean.grid(row=13,column=0, padx=2, pady=2,columnspan=3)
        btnSubmit_dsclean.grid(row=14,column=0, padx=2, pady=2)
        btnSubmit_gclean.grid(row=14,column=1, padx=2, pady=2)
        btnSubmit_pclean.grid(row=14,column=2, padx=2, pady=2)
                        
                 
        
        
        
        




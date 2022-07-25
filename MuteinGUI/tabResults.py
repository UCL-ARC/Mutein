
import tkinter as tk
from functools import partial
import tabConfig as rs
from io import StringIO as sio
import Analysis
import pandas as pd
import os
from os.path import exists

def binaryToDataFrame(binary):
    resA = str(binary,"utf-8")
    lines = resA.split("\n")
    res = ""
    ret_lines = []
    started = False
    for line in lines:
        if started:
            if "DATAFRAME_END" not in line:
                ret_lines.append(line.strip())
            else:
                started = False
        elif "DATAFRAME_START" in line:
            started = True
    
    if len(ret_lines) == 0:
        return resA
    # turn it into a dataframe
    line_format = [ret_lines[0].strip()]
    cols = ret_lines[0].split(",")    
    dic_data = {}    
    for i in range(len(cols)):
        dic_data[i] = []    
    total_len = len(cols)
    
    for rl in ret_lines[1:]:
        if len(rl) > 2:#a little arbitrary
            line_format.append(rl[0].strip())
            vals = rl.split(",")
            for i in range(len(vals)):
                dic_data[i].append(vals[i])
            if len(vals) != total_len:
                print(vals)
        

    df = pd.DataFrame.from_dict(dic_data)
    df.columns = cols                        
    return df

def show_background(txtBox,txtDataset,txtGene,txtPdb,mode):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    if mode=="GENE":    
        ret = rs.RunScript("GENE_BACK",dataset+":"+gene+":"+pdb)
    else:
        ret = rs.RunScript(mode,dataset+":"+gene+":"+pdb)
    df = binaryToDataFrame(ret)
    txtBox.delete(1.0,tk.END)
    txtBox.insert(tk.END, df)        
    txtBox.insert(tk.END, "\n")    
    txtBox.insert(tk.END, ret)                                    

def show_plots(txtBox,txtDataset,txtGene,txtPdb,mode):        
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    if mode == "PDB":
        pdb = txtPdb.get().strip()
    else:
        pdb = ""
    ret_back = rs.RunScript("PDB_BACK",dataset+":"+gene+":"+pdb)
    ret_var = rs.RunScript("PDB_BM",dataset+":"+gene+":"+pdb)
    df_back = binaryToDataFrame(ret_back)
    df_var = binaryToDataFrame(ret_var)
    ana = Analysis.Analysis(df_back,df_var,dataset+":"+gene+":"+pdb)
    
    if mode == "PDB":
        ana.createDdgBackRid()
        ana.createDdgBackFromTo()
        ana.histAllBackground()
        ana.show()
    else:
        ana.createPdbSummary()        
        ana.show()    

def save_dataframes(txtBox,txtPath, txtDataset,txtGene,txtPdb):
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()    
    pdb = ""   
    path = txtPath.get().strip()
    genes = []     
    if gene == "" or gene.lower()== "x":
        ret = rs.RunScript("DS_GENES",dataset+":"+gene+":"+pdb)
        df = binaryToDataFrame(ret)
        try:
            genes = df["gene"]
        except:
            pass
    else:
        genes.append(gene)
        
    for gene in genes:
        txtBox.delete(1.0,tk.END)
        print("Retrieving data for",gene)
                
        out_path = path + dataset.lower() + "_" + gene.lower() + ".csv"
        out_path_var = path + dataset.lower() + "_" + gene.lower() + "_var.csv"
        if not exists(out_path):
            ret_back = rs.RunScript("GENE_BACK",dataset+":"+gene+":"+pdb)
            df_back = binaryToDataFrame(ret_back)
            try:
                df_back.to_csv(out_path,index=False)
                print("Saved to",out_path)
                txtBox.insert(tk.END, out_path)        
                txtBox.insert(tk.END, df_back)        
            except:
                print("Error saving",out_path)
                print(df_back)
        else:
            print("File already exists",out_path)
                        
        if not exists(out_path_var):
            ret_var = rs.RunScript("GENE_VAR",dataset+":"+gene+":"+pdb)
            df_var = binaryToDataFrame(ret_var)
            try:
                df_var.to_csv(out_path_var,index=False)
                print("Saved to",out_path_var)
                txtBox.insert(tk.END, out_path_var)                
                txtBox.insert(tk.END, df_var)        
            except:
                print("Error saving",out_path_var)
                print(df_var)
        else:
            print("File already exists",out_path_var)

        
        
        
     
                
        
class tabResults:

    def __init__(self,parent):
        self.parent = parent
        self.txtBox = None
        self.txtJob = None
                                                
    def createTab(self,clr,clra):                
        ### layout ###
        
        # frame
        frameLHS = tk.Frame(self.parent, width=50,bg=clr)
        frameLHS.pack(padx=5, pady=15, side=tk.LEFT)    
        
        frame = tk.Frame(self.parent, width=50,bg=clr)
        frame.pack(padx=5, pady=15, side=tk.LEFT)    
        
        
        rhs = tk.Label(frameLHS,bg=clra,text="")
        rhs.pack(padx=5, pady=15, side=tk.TOP)
                
        # LEFT HAND SIDE ##
        lhs = tk.Label(frameLHS,bg=clra,text="")
        lhs.pack(padx=5, pady=15, side=tk.TOP)
        
        
        
        header = tk.Label(lhs,bg=clra,text="Explore Mutein Results",width=60)                        
        self.txtBox = tk.Text(lhs,height=40,width=100,borderwidth=5)                
        self.txtBox.insert('end', "")
        scrollb = tk.Scrollbar(lhs,command=self.txtBox.yview)        
        self.txtBox['yscrollcommand'] = scrollb.set

        header.grid(row=0,column=0, padx=2, pady=2,columnspan=3)
        self.txtBox.grid(row=1,column=0,rowspan=10)
        scrollb.grid(row=1, column=1, sticky='nsew',rowspan=10)

                
        # RIGHT HAND SIDE ##
        header2 = tk.Label(rhs,bg=clra,text="Results to choose",width=60)                        

                
        lblDataset = tk.Label(rhs,bg=clra,text="dataset",width=20)
        lblGene = tk.Label(rhs,bg=clra,text="gene",width=20)
        lblPdb = tk.Label(rhs,bg=clra,text="pdb",width=20)        
        text_dataset = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")        
        text_gene = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")        
        text_pdb = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")
        text_dataset.insert(0, "sheartop")
        text_gene.insert(0, "NOTCH1")
        text_pdb.insert(0, "1pb5")

                
        btnGene = tk.Button(rhs, text="Show Plots-Gene", command=partial(show_plots,self.txtBox,text_dataset,text_gene,text_pdb,"GENE"),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        btnPdb = tk.Button(rhs, text="Show Plots-Pdb", command=partial(show_plots,self.txtBox,text_dataset,text_gene,text_pdb,"PDB"),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        
        btnBG = tk.Button(rhs, text="View background ddg", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"PDB_BACK"),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        btnBM = tk.Button(rhs, text="View variant ddg", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"PDB_BM"),borderwidth=5, relief="groove",width=20,bg="goldenrod")

        # save the dataframes to the path given
        lblPath = tk.Label(rhs,bg=clra,text="local path",width=20)        
        text_path = tk.Entry(rhs,width=20,borderwidth=5, relief="sunken")                
        filepatha = os.path.dirname(os.path.realpath(__file__))
        filepath = filepatha.replace("\\","/")
        dir = os.path.dirname(filepath).split("/")[:-1]
        retpath = "/".join(dir) + "/Mutein/Results/"
        text_path.insert(0,retpath)
        btnPath = tk.Button(rhs, text="Save Dataset data", command=partial(save_dataframes,self.txtBox,text_path, text_dataset,text_gene,text_pdb),borderwidth=5, relief="groove",width=20,bg="coral")
                        
        header.grid(row=0,column=0, padx=2, pady=2,columnspan=3)        
        btnGene.grid(row=1,column=0, padx=2, pady=2)                
        btnPdb.grid(row=1,column=1, padx=2, pady=2)                
        btnBG.grid(row=2,column=0, padx=2, pady=2)                
        btnBM.grid(row=2,column=1, padx=2, pady=2)              
            
        lblDataset.grid(row=3,column=0, padx=2, pady=2)    
        lblGene.grid(row=3,column=1, padx=2, pady=2)
        lblPdb.grid(row=3,column=2, padx=2, pady=2)        
        text_dataset.grid(row=4,column=0, padx=2, pady=2)      
        text_gene.grid(row=4,column=1, padx=10, pady=2)
        text_pdb.grid(row=4,column=2, padx=2, pady=2)

        btnPath.grid(row=5,column=0, padx=2, pady=2)                  
        lblPath.grid(row=5,column=1, padx=2, pady=2)                  
        text_path.grid(row=5,column=2, padx=2, pady=2)                  
        
                        
                 
        
        
        
        




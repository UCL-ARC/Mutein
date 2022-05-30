
import tkinter as tk
from functools import partial
import tabConfig as rs
from io import StringIO as sio
import Analysis
import pandas as pd

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

def show_background(txtBox,txtDataset,txtGene,txtPdb,mode,tab,frame):    
    #try: 
    #    figure_canvas.get_tk_widget().pack_forget()
    #    figure_canvas.destroy()
    #    tab.destroy()
    #    tab = tk.Frame(frame, width=50)
    #    tab.pack(padx=5, pady=15, side=tk.TOP)    
    #except: 
    #    pass 
    dataset = txtDataset.get().strip()
    gene = txtGene.get().strip()
    pdb = txtPdb.get().strip()
    if mode=="GENE":    
        ret = rs.RunScript("PDB_BACK",dataset+":"+gene+":"+pdb)
    else:
        ret = rs.RunScript(mode,dataset+":"+gene+":"+pdb)
    df = binaryToDataFrame(ret)
    txtBox.delete(1.0,tk.END)
    txtBox.insert(tk.END, df)        
    txtBox.insert(tk.END, "\n")    
    txtBox.insert(tk.END, ret)                                    
    
    #ana = Analysis.Analysis(df)
    #if mode=="GENE":
    #    ana.createPdbSummary()
    #elif mode=="PDB_BACK":
    #    ana.createDdgBackResidue()
    
    
    

        
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

                
        btnCOV = tk.Button(rhs, text="View gene coverage", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"GENE",frame,frameLHS),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        btnBG = tk.Button(rhs, text="View background ddg", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"PDB_BACK",frame,frameLHS),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        btnBM = tk.Button(rhs, text="View variant bm ddg", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"PDB_BM",frame,frameLHS),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        btnPS = tk.Button(rhs, text="View variant ps ddg", command=partial(show_background,self.txtBox,text_dataset,text_gene,text_pdb,"PDB_PS",frame,frameLHS),borderwidth=5, relief="groove",width=20,bg="goldenrod")
        
        header.grid(row=0,column=0, padx=2, pady=2,columnspan=3)
        btnCOV.grid(row=1,column=0, padx=2, pady=2)                
        btnBG.grid(row=2,column=0, padx=2, pady=2)                
        btnBM.grid(row=2,column=1, padx=2, pady=2)                
        btnPS.grid(row=2,column=2, padx=2, pady=2)                
        lblDataset.grid(row=3,column=0, padx=2, pady=2)    
        lblGene.grid(row=3,column=1, padx=2, pady=2)
        lblPdb.grid(row=3,column=2, padx=2, pady=2)        
        text_dataset.grid(row=4,column=0, padx=2, pady=2)      
        text_gene.grid(row=4,column=1, padx=10, pady=2)
        text_pdb.grid(row=4,column=2, padx=2, pady=2)
        
                        
                 
        
        
        
        




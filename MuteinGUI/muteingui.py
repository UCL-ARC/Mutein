  
import tkinter as tk
from tkinter import ttk
from turtle import bgcolor

import tabMain
import tabConfig
import tabPipeline
import tabResults

class App():
    def __init__(self):
        super().__init__()

        clr_main = "dimgray"
        clr1,clr1a ="lightgreen","darkseagreen"
        clr2,clr2a="cornflowerblue","mediumslateblue"
        clr3,clr3a ="rosybrown","indianred"
        clr4,clr4a ="palegoldenrod","gold"


        rt = tk.Tk()
        rt.iconbitmap("C:/UCL/github/MuteinGUI/GUI/icon2.ico")
        rt.config(bg=clr2) 
        rt.title('QStat Viewer for Mutein')
        rt.geometry("1500x800")

        tab_parent = ttk.Notebook(rt)
        tab1 = tk.Frame(tab_parent,bg=clr1)
        tab2 = tk.Frame(tab_parent,bg=clr2)
        tab3 = tk.Frame(tab_parent,bg=clr3)
        tab4 = tk.Frame(tab_parent,bg=clr4)

        noteStyle = ttk.Style()
        noteStyle.theme_use('default')
        noteStyle.configure("TNotebook", background=clr_main, borderwidth=5)
        noteStyle.configure("TNotebook.Tab", background=clr_main, borderwidth=5)

        tab_parent.add(tab1,text="Config")
        tab_parent.add(tab2,text="Monitor")
        tab_parent.add(tab3,text="Pipeline")
        tab_parent.add(tab4,text="Results")
        tab_parent.pack(expand=1,fill="both")

        tC = tabConfig.tabConfig(tab1)
        tC.createTab(clr1,clr1a)
        tM = tabMain.tabMain(tab2)
        tM.createMainTab(clr2,clr2a)
        tP = tabPipeline.tabPipeline(tab3)
        tP.createTab(clr3,clr3a)
        tP = tabResults.tabResults(tab4)
        tP.createTab(clr4,clr4a)

        tk.mainloop()

    
if __name__ == '__main__':
    app = App()
    
    


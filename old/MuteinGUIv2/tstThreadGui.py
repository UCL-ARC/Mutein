"""
Code taken from: https://www.pythontutorial.net/tkinter/tkinter-thread/
"""
import tkinter as tk
from tkinter import ttk
from tkinter.messagebox import showerror
from threading import Thread
import SshComms as sc

class App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title('Webpage Download')
        self.geometry('1500x680')
        self.resizable(0, 0)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=5)

        self.create_view_frame()        
        self.create_header_frame()        
        self.create_body_frame()
        self.create_status_frame()
        self.create_footer_frame()
        self.create_tab_frame()

    def create_view_frame(self):
        self.view = ttk.Frame(self)
        self.view.grid(column=1, row=0, sticky=tk.NSEW, padx=10, pady=10)


    def create_tab_frame(self):
        clr_main = "dimgray"
        clr1,clr1a ="lightgreen","darkseagreen"
        clr2,clr2a="cornflowerblue","mediumslateblue"
        clr3,clr3a ="rosybrown","indianred"
        clr4,clr4a ="palegoldenrod","gold"
        self.tab_parent = ttk.Notebook(self)
        tab1 = tk.Frame(self.tab_parent,bg=clr1)
        tab2 = tk.Frame(self.tab_parent,bg=clr2)
        tab3 = tk.Frame(self.tab_parent,bg=clr3)
        tab4 = tk.Frame(self.tab_parent,bg=clr4)
        noteStyle = ttk.Style()
        noteStyle.theme_use('default')
        noteStyle.configure("TNotebook", background=clr_main, borderwidth=5)
        noteStyle.configure("TNotebook.Tab", background=clr_main, borderwidth=5)
        self.tab_parent.add(tab1,text="Config")
        self.tab_parent.add(tab2,text="Monitor")
        self.tab_parent.add(tab3,text="Pipeline")
        self.tab_parent.add(tab4,text="Results")
        #self.tab_parent.pack(expand=1,fill="both")

        self.tab_parent.grid(column=0, row=0, sticky=tk.NSEW, padx=10, pady=10)

        import tabConfig
        tC = tabConfig.tabConfig(self,tab1,self.html)
        tC.createTab(clr1,clr1a)
        import tabMonitor
        tM = tabMonitor.tabMonitor(self,tab2,tC,self.html)
        tM.createTab(clr2,clr2a)


    
    def create_header_frame(self):

        self.header = ttk.Frame(self.view)
        
        # attach the header frame
        self.header.grid(column=0, row=0, sticky=tk.NSEW, padx=10, pady=10)

    
    def monitor(self, thread):
        if thread.is_alive():
            # check the thread every 100ms
            self.after(100, lambda: self.monitor(thread))
        else:
            if thread.msg == "FLUSH":
                self.html.delete(1.0,tk.END)
                self.log.delete(1.0,tk.END)
            self.html.insert('end', thread.msg)                        
            self.log.insert('end', thread.error)
            #self.download_button['state'] = tk.NORMAL

    def create_body_frame(self):
        self.body = ttk.Frame(self.view)
        # text and scrollbar
        self.html = tk.Text(self.body, height=25,width=120)
        self.html.grid(column=0, row=1)

        scrollbarV = ttk.Scrollbar(self.body,
                                  orient='vertical',
                                  command=self.html.yview)

        scrollbarV.grid(column=1, row=1, sticky=tk.NS)
        self.html['yscrollcommand'] = scrollbarV.set

        # attach the body frame
        self.body.grid(column=0, row=1, sticky=tk.NSEW, padx=10, pady=10)
    
    def create_status_frame(self):
        self.status = ttk.Frame(self.view)
        # text and scrollbar
        self.log = tk.Text(self.status, height=5,width=120)
        self.log.grid(column=0, row=2)

        scrollbarL = ttk.Scrollbar(self.status,
                                  orient='vertical',
                                  command=self.log.yview)

        scrollbarL.grid(column=1, row=2, sticky=tk.NS)
        self.log['yscrollcommand'] = scrollbarL.set

        # attach the body frame
        self.status.grid(column=0, row=2, sticky=tk.NSEW, padx=10, pady=10)

    def create_footer_frame(self):
        self.footer = ttk.Frame(self.view)
        # configure the grid
        self.footer.columnconfigure(0, weight=1)
        # exit button
        self.exit_button = ttk.Button(self.footer,
                                      text='Exit',
                                      command=self.destroy)

        self.exit_button.grid(column=0, row=0, sticky=tk.E)

        # attach the footer frame
        self.footer.grid(column=0, row=3, sticky=tk.NSEW, padx=10, pady=10)


if __name__ == "__main__":
    app = App()
    app.mainloop()
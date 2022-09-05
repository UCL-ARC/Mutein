
import sys
import tkinter as tk
from functools import partial

# TODO Need to use background threads
# https://www.pythontutorial.net/tkinter/tkinter-thread/

# Modules only needed for this class
        
class tabMonitor:

    def __init__(self,app,parent,cfg, txtView):
        self.app = app
        self.parent = parent
        self.txtBox = txtView
        self.txtJob = None
        self.cfg = cfg
                                                
    def createTab(self,clr,clra):                
        ### layout ###        
        rhs = tk.Label(self.parent,bg=clra,text="")
        rhs.pack(padx=5, pady=15, side=tk.LEFT)

        rhs_up = tk.Label(rhs,bg=clra,text="")
        rhs_up.pack(padx=5, pady=15, side=tk.TOP)

        rhs_down = tk.Label(rhs,bg=clra,text="")
        rhs_down.pack(padx=5, pady=15, side=tk.TOP)
        
        header2 = tk.Label(rhs_up,bg="silver",text="Controls")
        header2.grid(row=0,column=0,columnspan=6,sticky = 'EW')
        
        self.txtJob = tk.Entry(rhs_up,width=20,borderwidth=5, relief="sunken")        
        self.txtJob.insert('end', "xyz")
        self.txtCancel = tk.Entry(rhs_up,width=20,borderwidth=5, relief="sunken")
        

        btnQ = tk.Button(rhs_up, text="Refresh QStat", command=partial(getQstat),borderwidth=5, relief="groove",width=20,bg="lime")
        btnC = tk.Button(rhs_up, text="Clean", command=partial(run_clean,self.txtBox),borderwidth=5, relief="groove",width=20,bg="darkorange")
        btnCA = tk.Button(rhs_up, text="Delete all log files", command=partial(run_clean_all,self.txtBox),borderwidth=5, relief="groove",width=20,bg="tomato")
        btnE = tk.Button(rhs_up, text="Rerun Errors", command=partial(run_rerun,self.txtBox),borderwidth=5, relief="groove",width=20, bg="red")
        btnL = tk.Button(rhs_up, text="View Log", command=partial(run_view,self.txtBox,self.txtJob),borderwidth=5, relief="groove",width=20,bg="lime")
        btnX = tk.Button(rhs_up, text="Cancel Jobs", command=partial(run_cancel,self.txtBox,self.txtCancel),borderwidth=5, relief="groove",width=20,bg="tomato")
        lblA = tk.Label(rhs_up,text="Formats accepted as patterns for the job delete:",width=60,bg=clra)
        lblB = tk.Label(rhs_up,text="Single job = 12345",width=60,bg=clra)
        lblC = tk.Label(rhs_up,text="Single task = 12345.7",width=60,bg=clra)
        lblD = tk.Label(rhs_up,text="From-to jobs = 12345-3456",width=60,bg=clra)
        lblE = tk.Label(rhs_up,text="From-to tasks = 12345.5-12345.9",width=60,bg=clra)
        lblF = tk.Label(rhs_up,text="!!! All !!! = *",width=60,bg="tomato")
        
        btnQ.grid(row=1,column=3, padx=2, pady=2)
        btnC.grid(row=2,column=3, padx=2, pady=2)
        btnCA.grid(row=2,column=4, padx=2, pady=2)
        btnE.grid(row=3,column=3, padx=2, pady=2)
        self.txtJob.grid(row=4,column=3, padx=2, pady=2)
        btnL.grid(row=4,column=4, padx=2, pady=2)
        self.txtCancel.grid(row=5,column=3, padx=2, pady=2)  
        btnX.grid(row=5,column=4, padx=2, pady=2)
        lblA.grid(row=6,column=3, padx=2, pady=2,columnspan=3)
        lblB.grid(row=7,column=3, padx=2, pady=2,columnspan=3)
        lblC.grid(row=8,column=3, padx=2, pady=2,columnspan=3)
        lblD.grid(row=9,column=3, padx=2, pady=2,columnspan=3)
        lblE.grid(row=10,column=3, padx=2, pady=2,columnspan=3)
        lblF.grid(row=11,column=3, padx=2, pady=2,columnspan=3)
        
                                                
        header2.grid(row=0,column=0,columnspan=6,sticky = 'EW')   

    def handle_ssh(self):   
        import SshComms as sc     
        un,pw,svr = self.cfg.getDetails()
        download_thread = sc.AsyncSsh()
        download_thread.setDetails(un,pw,svr)
        download_thread.start()
        self.app.monitor(download_thread)        


    def getQstat(self):
        self.cfg.runSshShell(["qstat"])    
    

    def cancelQJobs(self,job_task_pattern):  
        # possibilities to pass in
        
        # *
        # 1234-1239
        # 1234.2-1234.89
        # 1234.2
        # if tasks they cannot be 1234.2-1237.89, makes no sense
        
        #qdel JOB_ID -t 4-10  
        #qdel JOB_ID
        
        jobs_to_kill = []
        jobs = job_task_pattern.split("-")
        if job_task_pattern=="*":    
            jobs_to_kill.append('"*"')
        if len(jobs)==1:
            jobss = job_task_pattern.split(".")
            if len(jobss) == 1:
                jobs_to_kill.append(jobss[0])
            elif len(jobss) == 2:
                jobs_to_kill.append(jobss[0] + "-t " + str(jobss[1]))

        elif len(jobs)==2:
            from_job = jobs[0]
            to_job = jobs[1]
            jobssa = from_job.split(".")
            jobssb = to_job.split(".")
            if len(jobssa) == 2 and len(jobssb)==2:
                jobs_to_kill.append(jobssa[0] + "-t " + str(jobssa[1]) + "-" + str(jobssb[1]))
            elif len(jobssa) == 1 and len(jobssb)==1:
                for i in range(int(from_job),int(to_job)+1):
                    jobs_to_kill.append(str(i))

        results = []
        if len(jobs_to_kill)==0:
            return "You didn't enter a valid job-kill-pattern"
        
        for jobkill in jobs_to_kill:
            shell = self.getSshShell()
            try:    
                result = shell.run(["qdel",jobkill])
                results.append(result.output)
            except:
                error = sys.exc_info()[0]
                results.append(str(error) + " : " + jobkill)        
        
        return results

    def run_clean(self,txtBox):
            res = self.RunClean()
            txtBox.delete(1.0,tk.END)
            txtBox.insert('end', res)

    def run_clean_all(self,txtBox):
            res = self.RunCleanAll()
            txtBox.delete(1.0,tk.END)
            txtBox.insert('end', res)

    def run_view(self,txtBox,txtJob):
            pattern = txtJob.get().strip()
            res = self.cfg.RunView(pattern)
            txtBox.delete(1.0,tk.END)
            txtBox.insert('end', res)

    def run_rerun(self):
            res = self.RunRerun()
            self.txtBox.delete(1.0,tk.END)
            self.txtBox.insert('end', "rerun")

    def run_cancel(self,txtBox,txtCancel):
            pattern = txtCancel.get().strip()
            res = self.cancelQJobs(pattern)
            txtBox.delete(1.0,tk.END)
            txtBox.insert('end', res)
            


   
        
        



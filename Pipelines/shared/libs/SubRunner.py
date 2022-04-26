"""
RSA 21/4/22
------------------------
Class to manage the submission to of any out-of-process batch commands

"""
import os
import subprocess


class SubRunner:
    def __init__(self, exe, curr_dir,work_dir,script,ext,inputs,isarray):
        # need to add the array btach option here
        self.work_dir = work_dir
        self.args = []
        self.args.append(exe) #0 arg=executable        
        if not isarray:
            sh_script_name = "pipeline_array.sh"        
        else:
            sh_script_name = "pipeline_single.sh"                  
        py_script = curr_dir + work_dir + script + ".py"
        
        if exe == "bash":
            exe_script = curr_dir + sh_script_name
            print("### SubRunner: chmod +x",exe_script)
            os.system("chmod +x " + exe_script)
            self.args.append(exe_script) #1 arg=executable script     
        else:
            exe_script = curr_dir + work_dir + sh_script_name        
            self.args.append(py_script) #1 arg=executable script          
        
        self.args.append(inputs) #2 arg=inputs
        self.args.append(py_script) #3 arg=python script
        self.args.append(curr_dir+work_dir) #4 workspace
                        

    def run(self):
        print("### SubRunner.Run():",self.args)        
        print("### SubRunner.Run(): working dir", os.getcwd())        
        print("### SubRunner.Run() arguments",self.args)
        process = subprocess.Popen(
            args=self.args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

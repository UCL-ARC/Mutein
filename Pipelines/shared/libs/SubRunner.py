"""
RSA 21/4/22
------------------------
Class to manage the submission to of any out-of-process batch commands

"""
import os
import subprocess


class SubRunner:
    def __init__(self, exe, work_dir,script, inputs):
        self.work_dir = work_dir
        self.args = []
        if exe != "":  # it is a script if it is ""            
            self.args.append(exe)
        else:
            print("### SubRunner: chmod +x",script)
            os.system("chmod +x " + script)
        self.args.append(script)
        inputss = inputs.split(" ")
        for input in inputss:
            self.args.append(input)

    def run(self):
        print("### SubRunner.Run():",self.args)
        print("### SubRunner.Run(): changing directory to", self.work_dir)
        os.chdir(self.work_dir)
        print("### SubRunner.Run(): working dir", os.getcwd())        
        print("### SubRunner.Run() arguments",self.args)
        process = subprocess.Popen(
            args=self.args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

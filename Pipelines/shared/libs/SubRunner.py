"""
RSA 21/4/22
------------------------
Class to manage the submission to of any out-of-process batch commands

"""
import os
import subprocess


class SubRunner:
    def __init__(self, exe, script, inputs):
        self.args = []
        if exe != "":  # it is a script if it is ""            
            self.args.append(exe)
        else:
            print("chmod:",script)
            os.system("## SubRunner: chmod +x " + script)
        self.args.append(script)
        inputss = inputs.split(" ")
        for input in inputss:
            self.args.append(input)

    def run(self):
        process = subprocess.Popen(
            args=self.args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

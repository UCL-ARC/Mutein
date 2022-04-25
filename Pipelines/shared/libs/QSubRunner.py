"""
RSA 21/4/22
------------------------
Class to manage the submission to qsub on myriad

"""
import os
import subprocess


class QSubRunner:
    def __init__(self, script,work_dir, dependency, time, array, homeuser, inputs, print_only):
        os.system("chmod +x " + script)
        self.print_only = print_only
        self.args = []
        self.args.append("qsub")
        if str(dependency) != "-1":
            self.args.append("-hold_jid")
            self.args.append(dependency)
        if int(array) > 0:
            self.args.append("-t")
            self.args.append("1-" + str(array))
        self.args.append("-l")  # $ -l h_rt=5:00:0
        self.args.append("h_rt=" + time)
        self.args.append("-wd")  # $ -wd /home/ucbtlcr/Scratch/workspace
        self.args.append("/home/" + homeuser + "/Scratch/workspace")
        self.args.append(script)
        self.args.append(work_dir)
        self.args.append(inputs)

    def run(self):
        if self.print_only:
            print("### QSub.Run print only",self.args)
        else:
            process = subprocess.Popen(
                args=self.args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            result = process.communicate()
            print(result)  # e.g. Your job 588483 ("foldx-aggregate") has been submitted
            results = result[0].split(" ")
            jobid = results[2]
            if "." in jobid:
                results = jobid.split(".")
                jobid = results[0]
                return jobid
            return -1

"""
RSA 21/4/22
------------------------
Class to manage the submission to qsub on myriad

"""
import os
import subprocess


class QSubRunner:
    def __init__(
        self,
        runid,
        qsubid,
        script,
        install_dir,
        data_dir,
        pipe_dir,
        dependency,
        time,
        array,
        homeuser,
        inputs,
        cores,
        print_only,
    ):
        inputs += "@install_dir=" + install_dir
        inputs += "@data_dir=" + data_dir
        isarray = int(array) > 0
        if isarray:
            sh_script_name = install_dir + "Pipelines/foldx/libs/pipeline_array.sh"
        else:
            sh_script_name = install_dir + "Pipelines/foldx/libs/pipeline_single.sh"
        os.system("chmod +x " + sh_script_name)
        self.runid = runid
        self.print_only = print_only
        self.args = []
        self.args.append("qsub")  # 0 executable
        if str(dependency) != "-1":
            self.args.append("-hold_jid")
            self.args.append(dependency)
        if int(array) > 0:
            self.args.append("-t")
            self.args.append("1-" + str(array))
        # ask for cores from the config
        self.args.append("-pe")
        self.args.append("smp")
        self.args.append(cores)
        if script == "cleanup":  # redirect the ogs as we don't care about them
            self.args.append("-o")
            self.args.append("clean_out.txt")
            self.args.append("-e")
            self.args.append("clean_error.txt")
        self.args.append("-N")  # $ -N foldx-posscan
        self.args.append(qsubid)  # $
        self.args.append("-l")  # $ -l h_rt=5:00:0
        self.args.append("h_rt=" + time)
        self.args.append("-wd")  # $ -wd /home/ucbtlcr/Scratch/workspace
        self.args.append("/home/" + homeuser + "/Scratch/workspace")
        self.args.append(sh_script_name)  # 1 executable script
        self.args.append(inputs)  # 2 inputs
        self.args.append(install_dir + pipe_dir + script + ".py")  # 3 pyscript
        self.args.append(install_dir + pipe_dir)  # 4 workspace

    def run(self):
        if self.print_only:
            print("### QSub.Run() print only", self.args)
            return self.runid
        else:
            # print("### QSub.Run(): working dir", os.getcwd())

            fullcall = ""
            for arg in self.args:
                fullcall += arg + " "
            self.args.append(fullcall)

            print("### QSub.Run()", self.args)
            process = subprocess.Popen(
                args=self.args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            result = process.communicate()
            print(result)  # e.g. Your job 588483 ("foldx-aggregate") has been submitted
            if "job rejected" in result:
                return "x"
            else:
                results = result[0].split(" ")
                if len(results) > 1:
                    jobid = results[2]
                    if "." in jobid:
                        results = jobid.split(".")
                        jobid = results[0]
                        return jobid
                    return jobid
                return "x"

"""
RSA 21/4/22
------------------------
Class to manage the submission to of any out-of-process batch commands

"""
import os
import subprocess


class SubRunner:
    def __init__(
        self, exe, install_dir, data_dir, pipe_dir, script, ext, inputs, isarray
    ):
        # need to add the array btach option here
        inputs += "@install_dir=" + install_dir
        inputs += "@data_dir=" + data_dir
        self.args = []
        self.args.append(exe)  # 0 arg=executable
        if not isarray:
            sh_script_name = install_dir + "Pipelines/foldx/libs/pipeline_array.sh"
        else:
            sh_script_name = install_dir + "Pipelines/foldx/libs/pipeline_single.sh"
        py_script = install_dir + pipe_dir + script + ".py"

        if exe == "bash":
            exe_script = install_dir + pipe_dir + sh_script_name
            print("### SubRunner: chmod +x", exe_script)
            os.system("chmod +x " + exe_script)
            self.args.append(exe_script)  # 1 arg=executable script
        else:
            exe_script = install_dir + pipe_dir + sh_script_name
            self.args.append(py_script)  # 1 arg=executable script

        # self.args.append(py_script) #3 arg=python script
        self.args.append(inputs)  # 2 arg=inputs
        # self.args.append(install_dir + pipe_dir) #4 workspace

    def run(self):
        print("### SubRunner.Run():", self.args)
        print("### SubRunner.Run(): working dir", os.getcwd())
        print("### SubRunner.Run() arguments", self.args)
        process = subprocess.Popen(
            args=self.args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        result = process.communicate()
        print(result)
        return True

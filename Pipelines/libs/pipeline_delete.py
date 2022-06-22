"""
-----------------------------
RSA 16/05/22
-----------------------------

This aggregates the outputs from split positionscans into 1 file
-----------------------------
N.b this file may be run on the myriad clusters or on a local machine
-----------------------------
"""
import os
import pwd
import pandas as pd
from shutil import copyfile
from os.path import exists
import subprocess
import pathlib

# import from the shared library in Mutein/Pipelines/shared/lib
import sys

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
retpath = "/".join(dirs) + "/libs"
sys.path.append(retpath)
import Arguments
from os import listdir
from os.path import isfile, join
import os


def run_pipeline(args):    
    print(args)
    ##############################################
    ret = "Actions:"
    argus = Arguments.Arguments(args)
    mode = argus.arg("MODE")
    homeuser = pwd.getpwuid(os.getuid())[0]
    scratch_dir = "/home/" + homeuser + "/Scratch/workspace/"
    print("scratch_dir", scratch_dir)

    onlyfiles = [f for f in listdir(scratch_dir) if isfile(join(scratch_dir, f))]
    if mode == "ALL":
        print("### Delete all log files ###")
        for file in onlyfiles:
            if exists(scratch_dir + file):
                os.remove(scratch_dir + file)

    elif mode == "CLEAN":
        print("### Delete unnecessary log files ###")
        file_numbers = {}
        for file in onlyfiles:
            name = file.split(".")[0]
            number = file.split(".")[1][1:]
            if len(file.split(".")) > 2:
                number += "." + file.split(".")[2]

            if number not in file_numbers:
                file_numbers[number] = name

        for number, name in file_numbers.items():
            error_file = scratch_dir + name + ".e" + str(number)
            out_file = scratch_dir + name + ".o" + str(number)
            if exists(out_file) and exists(error_file):
                with open(error_file) as fr:
                    lines_err = fr.readlines()
                with open(out_file) as fr:
                    lines_out = fr.readlines()
                if len(lines_err) == 0:
                    if len(lines_out) > 0:
                        if lines_out[-1].strip() == "MUTEIN SCRIPT ENDED":
                            os.remove(error_file)
                            os.remove(out_file)
                            print("...removing", number, name)
                            ret += ",Clean:" + str(number)
                else:
                    print("Errors", number, name)
            else:
                print("Missing files", number, name)

    elif mode == "RERUNERROR" or mode == "RERUNORPHANS":
        print("### Rerun failed or stalled jobs ###")
        error_only = mode=="RERUNERROR"
        file_numbers = {}
        for file in onlyfiles:
            name = file.split(".")[0]
            number = file.split(".")[1][1:]
            if len(file.split(".")) > 2:
                number += "." + file.split(".")[2]

            if number not in file_numbers:
                file_numbers[number] = name

        for number, name in file_numbers.items():
            shall_rerun = False
            error_file = scratch_dir + name + ".e" + str(number)
            out_file = scratch_dir + name + ".o" + str(number)
            if exists(error_file):
                with open(error_file) as fr:
                    lines_err = fr.readlines()
            if len(lines_err) > 0:
                shall_rerun = True
            elif not error_only and exists(out_file):
                time = pathlib.Path(out_file).stat().st_mtime
                now = pd.Timestamp.now()
                hourdiff = (now-time).thrs
                print(hourdiff, out_file)

            
            
            if shall_rerun and exists(out_file):
                with open(out_file) as fr:
                    lines_out = fr.readlines()
                if len(lines_out) > 1:
                    cmd = lines_out[1].strip()
                    if len(cmd) > 3:
                        print("Rerunning cmd:", cmd)
                        # but if it is an array job we only want to run the single task
                        args = cmd.split(" ")
                        if "pipeline_array" in cmd:
                            script = args[-4]
                            # cmd = cmd.replace("pipeline_array","pipeline_single")
                            script = script[:-8]
                            script += "single.sh"
                            args[-4] = script
                            task = number.split(".")[1]
                            args[-3] += "@task=" + task
                            args = args[2:]
                            args[0] = "qsub"
                            print(number, name, "rerunning as single task:", args)

                        fullcall = ""
                        for arg in args:
                            fullcall += arg + " "
                        args.append(
                            fullcall
                        )  # add the batch itself onto the arguments

                        process = subprocess.Popen(
                            args=args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                        )
                        result = process.communicate()
                        print(result)
                        os.remove(error_file)
                        os.remove(out_file)
                        print("...removing", number, name)
                    else:
                        print("Can't re-run as manually started:", out_file)
            else:
                print("Can't re-run error as no outfile:", error_file)
    elif mode == "RERUNALL":
        print("### Rerun everything ###")
        file_numbers = {}
        for file in onlyfiles:
            name = file.split(".")[0]
            number = file.split(".")[1][1:]
            if len(file.split(".")) > 2:
                number += "." + file.split(".")[2]

            if number not in file_numbers:
                file_numbers[number] = name

        for number, name in file_numbers.items():
            error_file = scratch_dir + name + ".e" + str(number)
            out_file = scratch_dir + name + ".o" + str(number)
            if exists(out_file):
                with open(out_file) as fr:
                    lines_out = fr.readlines()
                if len(lines_out) > 1:
                    cmd = lines_out[1].strip()
                    if len(cmd) > 3:
                        print("Rerunning cmd:", cmd)
                        # but if it is an array job we only want to run the single task
                        args = cmd.split(" ")
                        if "pipeline_array" in cmd:
                            script = args[-4]
                            # cmd = cmd.replace("pipeline_array","pipeline_single")
                            script = script[:-8]
                            script += "single.sh"
                            args[-4] = script
                            task = number.split(".")[1]
                            args[-3] += "@task=" + task
                            args = args[2:]
                            args[0] = "qsub"
                            print(number, name, "rerunning as single task:", args)

                        fullcall = ""
                        for arg in args:
                            fullcall += arg + " "
                        args.append(fullcall)  # add the batch itself onto the arguments

                        process = subprocess.Popen(
                            args=args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                        )
                        result = process.communicate()
                        print(result)
                        os.remove(out_file)
                        print("...removing", number, name)
                        if exists(error_file):
                            os.remove(error_file)
            else:
                print("Can't re-run error as no outfile:", error_file)

    print("### COMPLETED CLEAN UP job ###")
    print("MUTEIN SCRIPT ENDED")
    return ret


##########################################################################################
if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)


import os
import sys
"""
This file sets up all the imports necessary at this level and any level below
The python files in python directory are the code gatekeepers
"""

dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
python_path = "/".join(dirs) + "/python"
dirs = os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]
lib_path = "/".join(dirs) + "/libs"
sys.path.append(python_path)
sys.path.append(lib_path)
#print(python_path)
#print(lib_path)
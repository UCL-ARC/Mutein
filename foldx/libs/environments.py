"""
RSA 11.4.22
--------------------------------------
New users, add your environment details here 
1) Add the 5 elements in the environment dictionary
2) Add to the IF statement a way to know who you are

"""
environments = (
    {}
)  
# [env]=[
#   (foldx path), 
#   (python path), 
#   (whether to run hpc or python),
#   (install_directory),
#   (data_directory)
# ]
environments["myriad"] = [
    "foldx", 
    "python", 
    "hpc",
    "",
    ""
]  # Main myriad environment

environments["rachel"] = [
    "~/UCL/libs/foldx5/foldx",
    "/bin/python3",
    "python",
    "/home/rachel/UCL/github/Mutein/",
    "/home/rachel/UCL/github/MuteinData/"
]  # RA laptop

environments["CI"] = [
    "~/UCL/libs/foldx5/foldx",
    "/bin/python3",
    "python",
    "",
    ""
]  # continuous integration

import os


def getenvironment(user=""):
    """Automatically recognise the environment though to can be overridden by exploicitly passing it in
    Returns: environment, tuple(foldx path, python path, hpc or python)
    """
    if user == "":
        dir_path = os.path.dirname(os.path.realpath(__file__)) + "/"
        if "/rachel/" in dir_path:
            user = "rachel"
        # elif add your user env here to default
        else:
            user = "myriad"
    if user in environments:
        envs = {}
        envs["foldxe"] = environments[user][0]
        envs["pythonexe"] = environments[user][1]
        envs["env"] = environments[user][2]
        envs["user"] = user
        envs["install_dir"] = environments[user][3]
        envs["data_dir"] = environments[user][4]
        return envs

    return {}

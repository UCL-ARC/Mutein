"""
RSA 11.4.22
--------------------------------------
New users, add your environment details here to and check them in to 
avoid having to mess about locally all the time
For example, when I run locally, I want to use python instead of hpc, 
and my exe path to foldx and python is different (environment variable incompetence)
"""
environments = (
    {}
)  # [env]=[triple of: (foldx path), (python path), (whether to run hpc or python or just log inputs)]
environments["myriad"] = ["foldx", "python", "hpc"]  # Main myriad environment
environments["myriad_tst"] = ["foldx", "python", "python"]  # run python tets on myriad
environments["rachel"] = [
    "~/UCL/libs/foldx5/foldx",
    "/bin/python3",
    "python",
]  # RA laptop
environments["inputs_hpc"] = [
    "foldx",
    "python",
    "inputs_hpc",
]  # just prints out what it would run
environments["inputs_python"] = [
    "foldx",
    "python",
    "inputs_python",
]  # just prints out what it would run
environments["python"] = [
    "foldx",
    "python",
    "python",
]  # just prints out what it would run
environments["CI"] = [
    "~/UCL/libs/foldx5/foldx",
    "/bin/python3",
    "python",
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
        return envs

    return {}

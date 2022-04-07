import sys, os


def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[
        : -1 * len("tests")
    ]  # add the path for the abocve directory
    sys.path.append(above_path)


ci_array_01 = [
    [
        "qsub",
        "-l",
        "h_rt=3:00:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx01_repair.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "1",
        "-l",
        "h_rt=0:05:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx02_makeparams.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "2",
        "-t",
        "1-50",
        "-l",
        "h_rt=6:00:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx03_posscan.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "3",
        "-l",
        "h_rt=0:05:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx04_aggddg.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "1",
        "-l",
        "h_rt=0:05:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx05_vparams.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "5",
        "-t",
        "1-63",
        "-l",
        "h_rt=1:00:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx06_build.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
    [
        "qsub",
        "-hold_jid",
        "6",
        "-l",
        "h_rt=0:05:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx07_vaggddg.sh",
        "6vxx",
        "6vxx_50",
        "50",
        ".",
        "Alpha",
        "6vxx_alpha_vars",
        "5",
    ],
]


def test_foldx_pipeline00_inputsA():
    addpath()
    import foldx00_pipeline as p00

    args = ["", "user=inputs_hpc", "pdb=6vxx"]
    ret_arr = p00.run_pipeline00(args)
    if ret_arr != ci_array_01:
        raise Exception("Fails CI 01")


ci_array_02 = [
    [
        "qsub",
        "-l",
        "h_rt=0:05:00",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx01_repair.sh",
        "1tst",
        "1tst_2",
        "2",
        ".",
        "Alpha",
        "1tst_vars",
        "2",
    ],
    [
        "qsub",
        "-hold_jid",
        "1",
        "-l",
        "h_rt=0:05:0",
        "-wd",
        "/home/rachel/Scratch/workspace",
        "foldx02_makeparams.sh",
        "1tst",
        "1tst_2",
        "2",
        ".",
        "Alpha",
        "1tst_vars",
        "2",
    ],
]


def test_foldx_pipeline00_inputsB():
    addpath()
    import foldx00_pipeline as p00

    args = [
        "",
        "user=inputs_hpc",
        "pdb=1tst",
        "jobs=12",
        "id=1@time=0:05:00",
        "id=3@time=2:00:00",
        "id=3@array=200",
    ]
    ret_arr = p00.run_pipeline00(args)
    if ret_arr != ci_array_02:
        raise Exception("Fails CI 02")


# test_foldx_pipeline00_inputsA()
test_foldx_pipeline00_inputsB()

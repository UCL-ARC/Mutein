"""
RSA 7/4/22
This runs the pipeline on a small dataset
"""
import sys, os


def addpath():
    above_path = os.path.dirname(os.path.realpath(__file__))[
        :-5
    ]  # TODO I will change this to use path helpers nect checkin
    sys.path.append(above_path)


def test_foldx_pipeline00_1tstAll():
    addpath()
    import foldx00_pipeline as p00

    args = ["", "user=CI", "jobs=1234567", "pdb=1tst", "repairs=3"]
    p00.run_pipeline00(args)


test_foldx_pipeline00_1tstAll()

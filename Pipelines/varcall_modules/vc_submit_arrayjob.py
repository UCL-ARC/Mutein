#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

import varcall as vc

def submit_arrayjob():
    '''
    get command from last line of file, insert an up to date timestamp
    submit the job to qsub
    '''
    args = parse_args()

    #skip to last line of file
    f = open(args.file)
    for line in f: pass
    f.close()

    #update the timestam[]
    cmd = line.strip().format(jobname_uid=vc.unique_id())

    #sanity check (qsub only for now)
    assert cmd.startswith("qsub")

    #issue the command, ensure no immediate error encountered
    #print(cmd)
    subprocess.run(cmd.split(),check=True)

def parse_args():
    'get the filename of the array job specification file'
    parser = argparse.ArgumentParser(description="Submit a qsub array job defined by an arrayjob specification file previously made by a generator")
    parser.add_argument('--file',   type=str, required=True, help='Path to the array job specification file')
    return parser.parse_args()

if __name__ == "__main__":
    globals()["submit_arrayjob"]()

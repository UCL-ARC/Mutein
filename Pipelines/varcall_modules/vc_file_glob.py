#!/usr/bin/env python

'''
generate tasklist file containing
'''

import glob
import os
import re
import sys
import json
import argparse
import varcall as vc

def generate_array_job():
    '''
    find the data files to send to fastqc by crawling through the "data" folder
    write the array job specification to a file one line per array task
    with the parameter names as a head and the qsub command as the footer
    '''

    args = parse_args()

    #setup new file with the per-task variable names listed in the header as column names
    t = vc.TaskManifest(args.taskfile,per_task=['accession'])

    ##crawl the data folder finding data files by globbing
    for accession in vc.glob_dirs(args.datadir,depth=3,**args.fileglob):
        t.add_task({"accession":accession})

    #write file footer containing the qsub command required to launch this array job
    t.close()

def parse_args():
    '''
    parse command line arguments and optional json file containing more arguments
    priority is command line options > last conf file... > first conf file > defaults
    '''

    parser = argparse.ArgumentParser(description='Generate qsub array job specification file for FastQC')
    parser.add_argument('--datadir',  type=str, help='folder containing datasets for FastQC to work on')
    parser.add_argument('--fileglob', type=json.loads, default={}, help='json string defining parameters for data directory globbing')
    parser.add_argument('--taskfile',   type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    vc.add_conf_arg(parser)
    args = vc.parse_and_load_conf(parser)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

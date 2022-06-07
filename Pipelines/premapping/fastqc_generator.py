#!/usr/bin/env python

'''
generate array job spec file for running FastQC on all selected samples
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
    f = vc.ArrayJob(args.out,fixed={"cores":args.cores},per_task=['accession'])

    ##crawl the data folder finding data files by globbing
    for accession in vc.glob_dirs(args.datadir,depth=3,**args.fileglob):
        f.write_task({"accession":accession})

    #write file footer containing the qsub command required to launch this array job
    f.write_qsub('fastqc',args)
    f.close()

def parse_args():
    '''
    parse command line arguments and optional json file containing more arguments
    priority is command line options > last conf file... > first conf file > defaults
    '''

    parser = argparse.ArgumentParser(description='Generate qsub array job specification file for FastQC')
    parser.add_argument('--datadir',   type=str, help='folder containing datasets for FastQC to work on')
    parser.add_argument('--fileglob',  type=json.loads, default={}, help='json string defining parameters for data directory globbing')
    parser.add_argument('--output',    type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    vc.add_standard_args(parser) #add grid engine related arguments
    args = vc.parse_and_load_conf(parser) #add any parameters from conf file(s)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

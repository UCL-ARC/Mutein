#!/usr/bin/env python

'''
generate array job spec file for running FastQC on all selected samples
'''

import glob
import os
import re
import sys
import argparse
import varcall as vc

def generate_array_job():
    '''
    write the array job specification to a file one line per array task
    with the parameter names as a head and the qsub command as the footer
    '''

    args = parse_args()

    #setup new file with the per-task variable names listed in the header as column names
    f = vc.ArrayJob(args.out,'threads accession')

    ##crawl the data folder finding data files by globbing
    for dataset in vc.glob_dirs(args.data,**args.dataset):
        for subset in vc.glob_dirs(dataset,**args.subset):
            for accession in vc.glob_dirs(subset,**args.accession):
                #write the next task to file as a new line
                f.write_task({"threads":args.threads,"accession":accession})

    #write file footer containing the qsub command required to launch this array job
    jobname = 'fastqc'
    f.write_qsub(jobname,args)
    f.close()

def parse_args():
    parser = argparse.ArgumentParser(description='Generate qsub array job specification file for FastQC')
    parser.add_argument('--data',     type=str, required=True, help='folder containing datasets for FastQC to work on')
    parser.add_argument('--out',      type=str, default='-',   help='file path to output array job spec to, - for stdout')
    parser.add_argument('--threads',  type=int, required=True, help='threads per FastQC job')
    parser.add_argument('--jsonfile', type=str,                help='json file defining optional globbing parameters for "dataset", "subset" and "accession"')
    vc.add_standard_args(parser)
    args =  parser.parse_args()

    #provide empty stubs for the optional globbing settings
    args.dataset = {}
    args.subset = {}
    args.accession = {}

    #load additional settings from JSON file
    if args.jsonfile:
        json_args = vc.load_json(args.jsonfile)
        if "dataset" in json_args: args.dataset = json_args["dataset"]
        if "subset" in json_args: args.subset = json_args["subset"]
        if "accession" in json_args: args.accession = json_args["accession"]  

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

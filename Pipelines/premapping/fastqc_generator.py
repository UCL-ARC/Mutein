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
    find the data files to send to fastqc by crawling through the "data" folder
    write the array job specification to a file one line per array task
    with the parameter names as a head and the qsub command as the footer
    '''

    args = parse_args()

    #setup new file with the per-task variable names listed in the header as column names
    f = vc.ArrayJob(args.out,'cores accession')

    ##crawl the data folder finding data files by globbing
    for dataset in vc.glob_dirs(args.data,**args.dataset):
        for subset in vc.glob_dirs(dataset,**args.subset):
            for accession in vc.glob_dirs(subset,**args.accession):
                f.write_task({"cores":args.cores,"accession":accession})

    #write file footer containing the qsub command required to launch this array job
    jobname = 'fastqc'
    f.write_qsub(jobname,args)
    f.close()

def parse_args():
    '''
    parse command line arguments
    '''
    parser = argparse.ArgumentParser(description='Generate qsub array job specification file for FastQC')
    parser.add_argument('--data',     type=str, help='folder containing datasets for FastQC to work on')
    parser.add_argument('--out',      type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    parser.add_argument('--jsonfile', type=str, help='json file defining optional globbing parameters for "dataset", "subset" and "accession"')
    parser.add_argument('--jsonstr',  type=str, help='json string defining optional globbing parameters for "dataset", "subset" and "accession"')
    vc.add_standard_args(parser) #add grid engine related arguments
    args =  parser.parse_args()

    #add default empty globbing filters
    args.dataset = {}
    args.subset = {}
    args.accession = {}

    #add optional arguments from JSON file and or string
    vc.add_json_arguments(args)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

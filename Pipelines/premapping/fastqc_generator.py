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
    parse command line arguments and optional json file containing more arguments
    priority is command line options > last json file... > first json file > defaults
    '''

    parser = argparse.ArgumentParser(description='Generate qsub array job specification file for FastQC')
    parser.add_argument('--data',      type=str, help='folder containing datasets for FastQC to work on')
    parser.add_argument('--out',       type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    parser.add_argument('--conf',      action='append',         help='json file(s) defining additional parameters')
    parser.add_argument('--dataset',   type=json.loads, default={}, help='json string defining parameters for "dataset" directory globbing')
    parser.add_argument('--subset',    type=json.loads, default={}, help='json string defining parameters for "subset" directory globbing')
    parser.add_argument('--accession', type=json.loads, default={}, help='json string defining parameters for "accession" directory globbing')
    vc.add_standard_args(parser) #add grid engine related arguments

    #parser values including conf file(s)
    args = vc.parse_and_load_conf(parser)

    print(args)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

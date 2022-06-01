#!/usr/bin/env python

'''
generate fastqc array job
'''

import glob
import os
import re
import sys
import argparse
import varcall as vc

def generate_array_job():
    args = parse_args()

    #setup new file with variable names as header
    f = vc.ArrayJob(args.output,'threads accession')

    threads=4 #this parameter stays the same across all tasks
    for dataset in vc.glob_dirs(args.datadir,exclude=['test$']):
        for subset in vc.glob_dirs(dataset,exclude=['gatk_db$','metadata$']):
            for accession in vc.glob_dirs(subset,exclude=[]):
                f.write_task(locals())#write "threads" and "accession" to next line of file

    jobname = 'fastqc'
    f.write_qsub(jobname,args)
    f.close()

def parse_args():
    parser = argparse.ArgumentParser(description='Generate qsub array job for FastQC of data')
    parser.add_argument('-d','--datadir', type=str, required=True, help='folder containing datasets as subfolders')
    parser.add_argument('-o','--output',  type=str, required=True, help='file to output array job spec to')
    vc.add_qsub_resource_params(parser)
    return parser.parse_args()

if __name__ == "__main__":
    globals()["generate_array_job"]()

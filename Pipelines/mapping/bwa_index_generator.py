#!/usr/bin/env python

'''
generate array job spec file preparing a reference genome for bwa and gatk
'''

import os
import sys
import argparse
import varcall as vc

def generate_array_job():
    '''
    this only deals with one reference genome sequence file at a time
    since we expect to use only one
    but you can change the path to point to any reference file
    '''

    args = parse_args()

    assert args.reference.endswith('.gz')

    f = vc.ArrayJob(args.out,'reference')
    f.write_task({"reference":args.reference})

    #write file footer containing the qsub command required to launch this array job
    jobname = 'bwa-index'
    f.write_qsub(jobname,args)
    f.close()

def parse_args():
    '''
    parse command line arguments
    '''

    parser = argparse.ArgumentParser(description='Generate qsub job specification file for indexing the reference sequence')
    parser.add_argument('--reference', type=str, help='path to the reference genome file which should be gzipped')
    parser.add_argument('--out',       type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    parser.add_argument('--jsonfile',  type=str, help='json file defining optional globbing parameters for "dataset", "subset" and "accession"')
    parser.add_argument('--jsonstr',   type=str, help='json string defining optional globbing parameters for "dataset", "subset" and "accession"')
    vc.add_standard_args(parser) #add grid engine related arguments
    args =  parser.parse_args()

    #add optional arguments from JSON file and or string
    vc.add_json_arguments(args)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

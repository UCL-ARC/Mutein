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

    f = vc.ArrayJob(args.out,fixed={'reference':args.reference})
    f.write_qsub('bwa-index',args)
    f.close()

def parse_args():
    '''
    parse command line arguments
    '''

    parser = argparse.ArgumentParser(description='Generate qsub job specification file for indexing the reference sequence')
    parser.add_argument('--reference', type=str, help='path to the reference genome file which should be gzipped')
    parser.add_argument('--output',    type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    vc.add_standard_args(parser) #add grid engine related arguments
    args = vc.parse_and_load_conf(parser)

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

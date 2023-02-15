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

    t = vc.TaskManifest(args.taskfile,per_task=['accession','read_group','read1','read2','bam'])
    f = vc.ArrayJob(args.out,fixed={'reference':args.reference})
    f.write_qsub('bwa-index',args)
    f.close()

def parse_args():
    '''
    parse command line arguments
    '''

    parser = argparse.ArgumentParser(description='Generate qsub job specification file for indexing the reference sequence')
    #parser.add_argument('--reference', type=str, help='path to the reference genome file which should be gzipped')
    parser.add_argument('--taskfile',    type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    vc.add_conf_arg(parser)
    args = vc.parse_and_load_conf(parser) 

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

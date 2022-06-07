#!/usr/bin/env python

'''
generate array job spec file for read mapping with bwa mem
'''

import os
import sys
import json
import argparse
import varcall as vc

def generate_array_job():
    '''
    discover fastq files to map against reference genome
    using bwa mem then to be sorted with samtools
    '''

    args = parse_args()

    #cores is the total cores, we give cores-2 to bwa mem, the remainder for samtools
    threads = int(args.cores) - 2
    assert threads > 0, "ensure at least three cores to account for samtools and bwa"

    f = vc.ArrayJob(args.out,fixed={'reference':args.reference,
                                    'cores':args.cores,
                                    'threads':threads},
                             per_task=['accession','read_group','read1','read2','bam'])

    ##crawl the data folder finding data files by globbing
    for accession in vc.glob_dirs(args.data,depth=3,**args.fileglob):
        read1 = vc.single_file(accession,remove_base=True,include=['_1.fastq.gz$'])
        read2 = vc.single_file(accession,remove_base=True,include=['_2.fastq.gz$'])
        read_group = '_'.join(accession.split('/')[1:])
        accession_id = accession.split('/')[-1]
        bam = accession_id+'_aln_sort.bam'
        f.write_task({"accession":accession,"read_group":read_group,
                        "read1":read1,"read2":read2,"bam":bam})

    #write file footer containing the qsub command required to launch this array job
    f.write_qsub('bwa-mem',args)
    f.close()

def parse_args():
    '''
    parse command line arguments
    '''

    parser = argparse.ArgumentParser(description='Generate qsub job specification file for bwa mem mapping reads to reference sequence and sorting them into a BAM file')
    parser.add_argument('--datadir',  type=str, help='folder to glob for read datasets')
    parser.add_argument('--fileglob', type=json.loads, default={}, help='json string defining parameters for data directory globbing')
    parser.add_argument('--output',   type=str, default='-',   help='file path to output array job spec to, - for stdout for debugging purposes')
    vc.add_standard_args(parser) #add grid engine related arguments
    args = vc.parse_and_load_conf(parser) 

    return args

if __name__ == "__main__":
    globals()["generate_array_job"]()

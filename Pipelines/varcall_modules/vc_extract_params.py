#!/usr/bin/env python

import os
import argparse
import json

def verify_params(header,params):
    '''
    verify that header and params contain same list of parameter names
    order need not be the same
    '''

    params_list = [item.strip() for item in params.strip().split()]
    params_list.sort()

    #do not alter the order of the actual column headers!
    header_list = [item for item in header]
    header_list.sort()

    assert header_list == params_list,"expected and provided parameters don't match"

def extract_params():
    args = parse_args()
    taskid = int(os.environ['SGE_TASK_ID'])

    f = open(args.tasklist)

    #get fixed parameter names and values
    fixed = json.loads(f.readline().strip())

    #get header (per-task parameter names)
    header = [item.strip() for item in f.readline().strip().split(',')]

    #parameter lists must match though need not be in same order
    verify_params(header+list(fixed.keys()),args.params)

    #skip to revelant line
    for i in range(taskid): line = f.readline()
    f.close()

    tokens = line.strip().split(',')
    assert len(header) == len(tokens)

    output = []

    #fixed parameter values
    for key,value in fixed.items():
        output.append('{key}="{value}"'.format(key=key,value=value))

    #per-task parameter values
    for i,key in enumerate(header):
        value = tokens[i]
        output.append('{key}="{value}"'.format(key=key,value=value))

    #set bash variables, assumes we are called using "$(vc_extract_params ...)"
    print('export ' + ' '.join(output))

def parse_args():
    parser = argparse.ArgumentParser(description="Extract the task's parameter values from the file, verify requested parameters match those in the file's header")
    parser.add_argument('--params',   type=str, required=True, help='space separated list of parameter names eg "threads accession"')
    parser.add_argument('--tasklist', type=str, required=True, help='tasklist file made by generator script')
    return parser.parse_args()

if __name__ == "__main__":
    globals()["extract_params"]()

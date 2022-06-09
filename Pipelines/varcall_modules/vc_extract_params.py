#!/usr/bin/env python

import os
import argparse
import json

def verify_params(header,params):
    '''
    verify that header and params contain same list of parameter names
    order need not be the same
    '''

    #params_list = [item.strip() for item in params.strip().split()]
    params_list = params
    params_list.sort()

    #do not alter the order of the actual column headers!
    header_list = [item for item in header]
    header_list.sort()

    assert header_list == params_list,"expected and provided parameters don't match"

def extract_params():
    args = parse_args()
    taskid = int(os.environ['SGE_TASK_ID'])

    if os.environ['vc_taskfile'] == 'no-task-file': return

    f = open(os.environ['vc_taskfile'])

    #get header (per-task parameter names)
    header = json.loads(f.readline().strip())

    #parameter lists match though need not be in same order
    verify_params(header,args.params)

    output = []

    #find the per-task parameters if present
    if len(header) > 0:
        #skip to relevant line
        count = 0
        for line in f:
            count += 1
            if count == taskid:
                break
        assert count == taskid, "insufficient lines in the tasklist file"
        tokens = line.strip().split(',')
        assert len(header) == len(tokens)

        #per-task parameter values
        for i,key in enumerate(header):
            value = tokens[i]
            output.append('{key}={value}'.format(key=key,value=value))

    f.close()

    #set bash variables, assumes we are called using "$(vc_extract_params ...)"
    if len(output) > 0:
        print('export ' + ' '.join(output))

def parse_args():
    parser = argparse.ArgumentParser(description="Extract the task's parameter values from the file, verify requested parameters match those in the file's header")
    parser.add_argument('params',  nargs='*', help='parameter names eg vc_extract_params threads accession')
    #parser.add_argument('--tasklist', type=str, default=None, help='tasklist file made by generator script')
    return parser.parse_args()

if __name__ == "__main__":
    globals()["extract_params"]()

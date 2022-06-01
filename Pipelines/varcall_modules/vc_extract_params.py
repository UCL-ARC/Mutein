#!/usr/bin/env python

import sys
import os

def extract_params():
    joblist = sys.argv[1]
    taskid = int(os.environ['SGE_TASK_ID'])

    f = open(joblist)
    header = f.readline().strip().split(',')
    for i in range(taskid): line = f.readline()
    f.close()

    tokens = line.strip().split(',')
    assert len(header) == len(tokens)

    output = []
    for i,key in enumerate(header):
        value = tokens[i]
        output.append('{key}={value}'.format(**locals()))
    print('export ' + ' '.join(output))

if __name__ == "__main__":
    globals()["extract_params"]()

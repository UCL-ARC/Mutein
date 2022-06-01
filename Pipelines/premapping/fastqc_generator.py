#!/usr/bin/env python

import glob
import os

'''
source ~/.mutein_settings

if [ $# -eq 0 ]
then
    JOBLIST_BASE=fastqc_joblist
else
    JOBLIST_BASE=$1
fi

rm -f ${JOBLIST_BASE}
'''

def list_dirs(path):
    '''
    return non-recursive list of subfolders of <path>
    '''
    next(os.walk(path))[1]

def generate_lines(filename,sep=None,comment=None,skip=0,strip=True):
    '''
    generator function which returns each line as a string
    return each line as a list of tokens if sep if not None
    ignore lines starting with the comment string
    use skip to skip any a fixed number of non-commented header lines
    strip surrounding whitespace unless strip is False
    '''

    f = open(filename)
    skipped = 0
    for line in f:
        #ignore commented lines
        if comment is not None and line.startswith(comment): continue

        #skip header lines
        if skipped < skip:
            skipped += 1
            continue

        #return tokens
        if sep is not None:
            if strip:
                yield line.strip().split(sep)
            else:
                yield line.split(sep)

        #return whole string
        if strip:
            yield line.strip()
        else:
            yield line

    f.close()

datadir = os.getcwd()

def recursive_filtered_glob(curr_path,curr_depth,ops):
    '''
    recursive glob
    '''

    for item in glob.iglob("**",root_dir=curr_path,recursive=False):
        new_path = os.path.join(curr_path,item)

        #if we need to recurse into subfolders
        if curr_depth < ops["max"]:
            if os.path.isdir(new_path):
                recursive_filtered_glob(new_path,curr_depth+1,min_depth)

        #if we need to yield items
        if curr_depth >= min_depth and curr_depth <= max_depth:
            if ops["dirs"] and os.path.isdir(new_path):
                yield new_path
            elif ops["files"] and os.path.isfile(new_path):
                yield new_path
            #if ops["syms"] and os.path.islink(new_path): yield new_path #islink and isdir/isfile can both be true

def glob_items(path,depth=0,
               include=[],exclude=[],include_files=[],exclude_files=[],
               dirs=True,files=True,symlinks=True):
    '''
    glob items under a given path with optional recursion
    optional include and exclude pattern strings
    optional include and exclude pattern files (containing patterns)
    include or exclude directories, files, symlinks
    filters block entering of folders as well as return items
    depth can be integer or tuple (min depth, max depth)
    '''

    #load filter files and add them to the include / exclude lists
    for fname in include_files: include += [line for line in generate_lines(fname)]
    for fname in exclude_files: exclude += [line for line in generate_lines(fname)]

    if type(depth) == int:
        min_depth = max_depth = depth
    else:
        min_depth,max_depth = depth

    curr_depth = 0
    ops = {"inc":include,
           "exc":exclude,
           "dirs":dirs,
           "files":files,
           "syms":symlinks,
           "min":min_depth,
           "max":max_depth,
          }

    for item in recursive_filtered_glob(path,curr_depth,ops):
        yield item



for dataset in glob_dirs(datadir):


'''
for DATASET in $(cat datasets/active_datasets)
do
    for SUBSET in $(cat datasets/${DATASET}/active_subsets)
    do
        for ACCESSION in $(cat datasets/${DATASET}/${SUBSET}/active_accessions)
        do
            echo "fastqc -t 4 datasets/${DATASET}/${SUBSET}/${ACCESSION}/*.fastq.gz" >> ${JOBLIST_BASE}
        done
    done
done

TOTAL_TASKS=$(cat ${JOBLIST_BASE} | wc --lines)

#for now just print qsub command to console for user to manually paste and run
#but we could call it directly, eg wrapped it in a screen session,
#use -sync y to catch errors etc
#eg screen -S fastqc -d -m qsub -sync y...
JOBNAME=fastqc-$(mutein_random_id)
echo qsub -N ${JOBNAME} -t 1-${TOTAL_TASKS} ${MUT_DIR}/Pipelines/premapping/fastqc_runner.sh ${JOBLIST_BASE}
echo capture_qacct.sh ${JOBNAME}
'''

#!/usr/bin/env python

'''
refactor of fastqc_generator.sh to use python globbing to determine datasets to process
rather than bash with simple list files

stuff not yet refactored into python:
source ~/.mutein_settings

if [ $# -eq 0 ]
then
    JOBLIST_BASE=fastqc_joblist
else
    JOBLIST_BASE=$1
fi

rm -f ${JOBLIST_BASE}
'''

import glob
import os
import re
import sys

def run_pipeline(args):
    datadir = args[1]

    for dataset in glob_dirs(datadir,depth=(0,2),include=['yokoyama2019/']):
        print(dataset)

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


def exclude_item(item,ops):
    '''
    check item against the exclusion and inclusion filters
    return True to exclude, False to include
    '''

    #apply exclusion filters
    for filt in ops["exc"]:
        if filt.search(item) is not None:
            #matches an exclusion filter, reject item
            return True

    #accept item if no inclusion filters provided
    if len(ops["inc"]) == 0: return False

    for filt in ops["inc"]:
        if filt.search(item) is not None:
            #matches an inclusion filter, accept item
            return False

    #didn't match an inclusion filter, reject item
    return True

def recursive_filtered_glob(curr_path,curr_depth,ops):
    '''
    recursive glob with filtering
    use non-recursive iglob call to iterate the current level only
    then selectively descend into sub folders by calling itself again
    idea: block enumeration of unwanted folders altogether
    to allow better scaling to millions of files and folders
    '''

    for item in glob.iglob(os.path.join(curr_path,"*"),recursive=False):
        #add trailing separator to allow regex matching of directories versus files
        if os.path.isdir(item): item += os.sep

        #descend into subfolders
        if os.path.isdir(item):
            if curr_depth < ops["max"]:
                for subitem in recursive_filtered_glob(item,curr_depth+1,ops):
                    yield subitem

        #apply filters
        if exclude_item(item,ops): continue

        #yield the current item
        if curr_depth >= ops["min"] and curr_depth <= ops["max"]:
            if ops["dirs"] and os.path.isdir(item):
                yield item
            elif ops["files"] and os.path.isfile(item):
                yield item

def glob_files(path,depth=0,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only files'
    for item in glob_items(path,depth,False,True,include,exclude,include_files,exclude_files):
        yield item

def glob_dirs(path,depth=0,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only dirs'
    for item in glob_items(path,depth,True,False,include,exclude,include_files,exclude_files):
        yield item

def glob_items(path,depth=0,dirs=True,files=True,
               include=[],exclude=[],include_files=[],exclude_files=[]):
    '''
    glob items under a given path with optional recursion
    optional include and exclude pattern strings
    optional include and exclude pattern files (containing patterns)
    include or exclude directories, files (symlinks treated as the thing they point to)
    filters block entering of folders as well as return of items
    depth can be integer (the min and max) or tuple (min depth, max depth)
    '''

    #load filter files and add them to the include / exclude lists
    for fname in include_files: include += [line for line in generate_lines(fname)]
    for fname in exclude_files: exclude += [line for line in generate_lines(fname)]

    #compile filters into regular expression objects
    include = [re.compile(filt) for filt in include]
    exclude = [re.compile(filt) for filt in exclude]

    if type(depth) == int:
        min_depth = max_depth = depth
    else:
        min_depth,max_depth = depth

    ops = {"inc":include,"exc":exclude,
           "dirs":dirs,"files":files,
           "min":min_depth,"max":max_depth}

    #print(ops)

    for item in recursive_filtered_glob(path,0,ops):
        yield item

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

if __name__ == "__main__":
    import sys

    globals()["run_pipeline"](sys.argv)

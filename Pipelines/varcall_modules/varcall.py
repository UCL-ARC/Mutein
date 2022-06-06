import glob
import os
import errno
import re
import sys
import random
import time
import json
from datetime import datetime

#ensure random has been seeded from system urandom or time (is this already done by the import??)
random.seed()

def unique_id():
    'generate unique id based on microsecond timestamp and pseudo random number'

    timestamp = datetime.fromtimestamp(time.time()).strftime("%Y%m%d-%H%M%S.%f")[:-3]
    randstamp = "-%04d"%random.randrange(1e4)

    return timestamp + randstamp

def find_file(filename):
    '''
    return path to the filename
    no change if path exists
    otherwise look in $MUT_CONFIG_DIR if defined
    otherwise raise exception
    '''

    if os.path.exists(filename): return filename

    if "MUT_CONFIG_DIR" in os.environ:
        filename = os.path.join(os.environ['MUT_CONFIG_DIR'],filename)

    if os.path.exists(filename): return filename

    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)

def add_json_arguments(args):
    '''
    load additional arguments from a json file and/or string to args
    --jsonfile <filename>
    --jsonstr <json string>
    do not overwrite existing arguments from the command line
    '''

    #load additional args from JSON file
    if args.jsonfile:
        filename = find_file(args.jsonfile)
        f = open(filename)
        for key,val in json.load(f).items():
            if not hasattr(args,key):
                setattr(args,key,val)
        f.close()

    #load additional args from JSON string
    if args.jsonstr:
        for key,val in json.loads(args.jsonstr).items():
            if not hasattr(args,key):
                setattr(args,key,val)

def add_standard_args(parser):
    '''
    add common job scheduler related parameters to the command line argument list
    intended to be generic, but based on grid engine
    '''

    parser.add_argument('-m','--mem',   type=str, default='4G',       help='memory per core <integer>[M,G,T]')
    parser.add_argument('-t','--time',  type=str, default='04:00:00', help='max wall time HH:MM:SS')
    parser.add_argument('-c','--cores', type=str, default='4',        help='cores')
    parser.add_argument('-p','--tmpfs', type=str, default='10G',      help='tmpfs space <integer>[M,G,T]')
    parser.add_argument('-l','--logs',  type=str, default='logs',     help='path to write stdout and stderr log files to, "none" to rely on default behaviour')
    return parser

class ArrayJob:
    def __init__(self,filename,header,sep=','):
        self.tasks = 0
        self.filename = filename

        #convert header into a list if it's given as a space delimited string
        if isinstance(header,str):
            self.header = header.strip().split()
        else:
            self.header = header

        self.sep = sep

        if self.filename == '-':
            self.f = sys.stdout
        else:
            self.f = open(self.filename,'w')

        self.f.write(self.sep.join(self.header)+'\n')
        self.format = self.sep.join(['{%s}'%x for x in self.header])

    def write_task(self,args):
        self.f.write(self.format.format(**args)+'\n')
        self.tasks += 1

    def append(self,line):
        self.f.write(line)
        if not line.endswith('\n'): self.f.write('\n')

    def write_qsub(self,name,args):
        jobname = name
        jobscript = os.path.join(os.environ['MUT_TEMPLATE_DIR'],name+'.qsub')

        cmd  = "qsub -cwd -V"

        if args.logs != "none":
            cmd += " -o {args.logs}/\$JOB_NAME.\$TASK_ID.o -e {args.logs}/\$JOB_NAME.\$TASK_ID.e"

        cmd += " -N {jobname}-{{jobname_uid}} -t 1-{self.tasks}"
        cmd += " -l h_rt={args.time} -l mem={args.mem} -l tmpfs={args.tmpfs} -pe smp {args.cores}"
        cmd += " {jobscript} {self.filename}"

        self.append(cmd.format(**locals()))

    def ntasks(self):
        return self.tasks

    def close(self):
        self.f.close()
        return self.tasks

    def __exit__(self,exception_type, exception_value, exception_traceback):
        self.f.close()

    def __enter__(self):
        return self

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

def glob_files(path,depth=1,remove_base=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only files'
    for item in glob_items(path,depth,False,True,remove_base,include,exclude,include_files,exclude_files):
        yield item

def glob_dirs(path,depth=1,remove_base=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only dirs'
    for item in glob_items(path,depth,True,False,remove_base,include,exclude,include_files,exclude_files):
        yield item

def glob_items(path,depth=1,dirs=True,files=True,remove_base=False,
               include=[],exclude=[],include_files=[],exclude_files=[]):
    '''
    glob items under a given path with optional recursion
    optional include and exclude pattern strings
    optional include and exclude pattern files (containing patterns)
    include or exclude directories, files (symlinks treated as the thing they point to)
    filters don't block recursion into folders only the return of items
    use depth to control excessive recursion or manually glob individual subfolders
    depth can be integer (the min and max) or tuple (min depth, max depth)
    '''

    #load filter from file and add them to the include / exclude lists
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

    for item in recursive_filtered_glob(path,1,ops):
        if remove_base == True:
            yield item[len(path):]
        else:
            yield item

def list_dirs(path):
    '''
    return non-recursive list of subfolders of <path>
    '''
    next(os.walk(path))[1]

def list_files(path):
    '''
    return non-recursive list of files in <path> folder
    '''
    next(os.walk(path))[2]

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

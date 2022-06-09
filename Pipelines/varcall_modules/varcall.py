import glob
import os
import errno
import re
import sys
import random
import time
import json
import argparse
from datetime import datetime

#ensure random has been seeded from system urandom or time (is this already done by the import??)
random.seed()

class JobManifest:
    '''
    contains the per-job parameters (ie those that stay the same in value between tasks)
    write a file containing the parameter names and values
    and the task manifest filename
    and a list of the jobs that can be launched
    '''
    def __init__(self,filename,jobname,args,task,fixed={}):
        'open the output file, write the parameter names and values'
        self.fixed = fixed
        self.filename = filename

        if self.filename == '-':
            self.f = sys.stdout
        else:
            self.filename += '.job'
            self.f = open(self.filename,'w')

        self.f.write(json.dumps(self.fixed)+'\n')

        task_count = task.n_tasks()
        taskfile = task.get_filename()
        jobscript = os.path.join(os.environ['MUT_TEMPLATE_DIR'],jobname+'.qsub')

        cmd   = "TASKFILE='{taskfile}' SELECTFILE='{{selectfile}}'"
        cmd  += " qsub -cwd -V"

        if args.logs != "none":
            cmd += " -o '{args.logs}/$JOB_NAME.$TASK_ID.o' -e '{args.logs}/$JOB_NAME.$TASK_ID.e'"

        cmd += " -N {jobname}-{{jobname_uid}} -t 1-{task_count}"
        cmd += " -l h_rt={args.time} -l mem={args.mem} -l tmpfs={args.tmpfs} -pe smp {args.cores}"
        cmd += " {jobscript}"

        self.job_list.append(cmd.format(**locals()))


class TaskManifest:
    '''
    contains the per-task parameters (ie those that differ in value between tasks)
    write a file containing a header of parameter names as a json list
    and one line per task of an array job containing the corresponding parameter values
    currently stored as csv style comma separated columns
    '''
    def __init__(self,filename,per_task=[]):
        'open the output file, write the header of parameter names'
        self.per_task = per_task
        self.format = ','.join(['{%s}'%x for x in self.per_task])
        self.task_count = 0
        self.filename = filename

        if self.filename == '-':
            self.f = sys.stdout
        else:
            self.f = open(self.filename,'w')

        self.f.write(json.dumps(self.per_task)+'\n')

    def add_task(self,kwargs):
        assert len(self.per_task) > 0, "no per-task parameters defined!"
        self.f.write(self.format.format(**kwargs)+'\n')
        self.task_count += 1

    def get_filename(self):
        return self.filename

    # def add_job(self,name,args):
    #     jobname = name

    #     #find job script template
    #     jobscript = os.path.join(os.environ['MUT_TEMPLATE_DIR'],name+'.qsub')

    #     cmd  = "qsub -cwd -V"

    #     if args.logs != "none":
    #         cmd += " -o '{args.logs}/$JOB_NAME.$TASK_ID.o' -e '{args.logs}/$JOB_NAME.$TASK_ID.e'"

    #     cmd += " -N {jobname}-{{jobname_uid}} -t 1-{self.tasks}"
    #     cmd += " -l h_rt={args.time} -l mem={args.mem} -l tmpfs={args.tmpfs} -pe smp {args.cores}"
    #     cmd += " {jobscript} {self.filename}"

    #     self.job_list.append(cmd.format(**locals()))

    def n_tasks(self):
        'how many tasks defined so far'
        return self.task_count

    # def njobs(self):
    #     'how many jobs defined so far'
    #     return len(self.job_list)

    def next_task_id(self):
        '''
        the number of the next task
        ie the task number that will be created by the next call to write_task
        '''
        return self.ntasks()+1

    # def write_manifest(self,filename):
    #     if filename == '-':
    #         self.f = sys.stdout
    #     else:
    #         sellf.f = open(filename,'w')

    #     self.f.write(json.dumps({"jobs":len(self.job_list),"tasks":len(self.task_list)})+'\n')
    #     for job in self.job_list: self.f.write(job + '\n')
    #     self.f.write(json.dumps(fixed)+'\n') #fixed parameter names and values
    #     self.f.write(json.dumps(self.per_task)+'\n') #names of per-task parameters
    #     self.f.close()

    def close(self):
        self.f.close()
        return self.task_count

    # def __exit__(self,exception_type, exception_value, exception_traceback):
    #     self.f.close()

    # def __enter__(self):
    #     return self

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

def parse_and_load_conf(parser):
    '''
    load additional arguments from any --conf option(s)
    store the values as the defaults in the parser
    ready to re-parse the command line options in order to give them priority
    '''

    #parse command line arguments
    args = parser.parse_args()
    if args.conf is None: return args
    if type(args.conf) == str: args.conf = [args.conf]

    #read conf files, store as new default values
    for jsonfile in args.conf:
        filename = find_file(jsonfile) # in case json is in the conf not the local folder
        f = open(filename)
        loaded = json.load(f)
        parser.set_defaults(**loaded)
        f.close()

    #re-read command line
    return parser.parse_args()

'''
class JSON(argparse.Action):
    'action to parse a json string, but found type=json.loads works better'
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super().__init__(option_strings, dest, **kwargs)
    def __call__(self, parser, namespace, values, option_string=None):
        print('%r %r %r' % (namespace, values, option_string))
        setattr(namespace, self.dest, json.loads(values))
'''

def add_conf_arg(parser):
    'add the conf file option'
    parser.add_argument('--conf',  action='append',         help='json file(s) defining additional parameters')
    return parser

def get_standard_args():
    'keep this in sync with the options in add_standard_args'
    return ['mem','time','cores','tmpfs','logs']

def add_standard_args(parser):
    '''
    add job scheduler related parameters to the command line argument list
    intended to be generic, but based on grid engine
    '''

    parser.add_argument('--mem',   type=str, default='4G',       help='memory per core <integer>[M,G,T]')
    parser.add_argument('--time',  type=str, default='01:00:00', help='max wall time HH:MM:SS')
    parser.add_argument('--cores', type=str, default='1',        help='cores')
    parser.add_argument('--tmpfs', type=str, default='10G',      help='tmpfs space <integer>[M,G,T]')
    parser.add_argument('--logs',  type=str, default='logs',     help='path to write stdout and stderr log files to, "none" to rely on default behaviour')
    return parser

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

def single_file(path,depth=1,remove_base=False,allow_spaces=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    return single_item(path,depth,False,True,remove_base,allow_spaces,include,exclude,include_files,exclude_files)

def single_dir(path,depth=1,remove_base=False,allow_spaces=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    return single_dir(path,depth,True,False,remove_base,allow_spaces,include,exclude,include_files,exclude_files)

def single_item(path,depth=1,dirs=True,files=True,remove_base=False,allow_spaces=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    '''
    convenience wrapper to glob_items returning a single item
    or raising exception for multiple hits
    '''

    hit = None
    for item in glob_items(path,depth,dirs,files,remove_base,allow_spaces,include,exclude,include_files,exclude_files):
        if hit != None:
            raise Exception('more than one matching item')
        else:
            hit = item

    if hit == None:
        raise Exception('no matching items found')

    return hit

def glob_files(path,depth=1,remove_base=False,allow_spaces=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only files'
    for item in glob_items(path,depth,False,True,remove_base,allow_spaces,include,exclude,include_files,exclude_files):
        yield item

def glob_dirs(path,depth=1,remove_base=False,allow_spaces=False,include=[],exclude=[],include_files=[],exclude_files=[]):
    'convenience wrapper to glob_items returning only dirs'
    for item in glob_items(path,depth,True,False,remove_base,allow_spaces,include,exclude,include_files,exclude_files):
        yield item

def glob_items(path,depth=1,dirs=True,files=True,remove_base=False,allow_spaces=False,
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
            item = item[len(path)+1:]

        item = os.path.normpath(item)
        
        if allow_spaces == False:
            assert not ' ' in item,'space(s) found in filename {item}'.format(item)

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

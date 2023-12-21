# YAMLmake Workflow Tool

This markdown documents YAMLmake, a new workflow tool developed during the Mutein project. This is similar to Snakemake in that a config file documents each step in the pipeline in terms of input and output files, the shell command(s) to be run and also allows local or HPC execution. However YAMLmake config files are written in pure YAML rather than a mixture of YAML and Snakemake exhanced Python. YAMLmake also does not create a dependency graph working backwards from a target output file, rather the actions are carried out sequentially in the order they appear in the file. This makes it much simpler to deal with lists of input files that are not known in advance, the action simply globs them before it runs, whereas Snakemake requires a special checkpoint concept to handle this common situation. YAMLmake allows easy changing between local and qsub execution, whereas Snakemake requires a separate directive for each rule to make it run remotely, and these can break when a module is imported. YAMLmake tries to specify all input and output files using special placeholders which insert filenames from predefined lists or from filesystem globs executed just prior to running the action. It therefore does not use embedded Python functions to generate the file lists.

## Key aims

- Produce a basic workflow engine that is easy to learn and not off-putting to new data scientists learning Linux and HPC on the job as part of their broader research
- No requirement to know any scripting language beyond the bash or script commands you would be running manually without a workflow engine
- Support reproducibility / portability through conda
- Support local and HPC execution through GridEngine
- Support array jobs on GridEngine
- Easily rerun only the failed jobs from an array
- Easily specify sets of input and output files for each action using placeholders rather than embedded Python functions
- Flexible, modular config files with comments to support well organised, maintainable pipelines

## Possible future features

- SLURM and Singularity / Apptainer support

## YAMLmake overview

The following outline of YAMLmake assumes some knowledge of bash to understand the details of the commands being run.

YAMLmake runs a YAML pipeline file containing a list of *actions* from top to bottom executing each one in the order encountered, taking account of *includes* and *modules*. Command line options allow jumping to run a specific action(s). A configuration state is also maintained as a hierachical tree of nested lists, dictionaries and strings which is modified everytime a *config* item is encountered in the pipeline. For example:

```
  - include: "base_config.yml"
  - action:
      name: "download_reference"
      exec: "local"
      input:
        url: "metadata/GRCh38_reference_url.txt"
      output:
        fastq: "reference/GRCh38.fasta.gz"
      shell: |
        wget $(cat {%url}) -O {%fastq}
```

Each action specifies one or more input and output files. If the output files are all present and more recent than any input file then the action is not executed. Likewise if any input file is missing the action will not run. The shell field of the action contains curly bracket placeholders where YAMLmake substitutes the input and output filenames before running the command either locally or through GridEngine. Provided the return code is zero and all the specified output files were created the action is considered to have succeeded. If the action is considered to have failed then any output files present will be removed, either by deletion or moving to a specified recycle bin folder. They can also be flagged as stale in-place by setting the mtime to a special value.

### include

The `include` item above causes the referred to file to be inserted into the pipeline as if it had been copy-pasted in at the location of the include. The included file can contain any YAMLmake configuration, including further nested includes and modules. When a pipeline first starts running the configuration state is set to build-in default values for the internal YAMLmake config fields. These then get overridden by any `config` items encountered in the main pipeline file or any includes.

### config

In the example above the file `base_config.yml` is included as the first item. This file would therefore contain the starting configuration options for the pipeline in the form of a `config` item, for example: 

```
- config:
    ym:
      prefix:               "my_project."
      stale_output_file:    "recycle"
      stale_output_dir:     "ignore"
      missing_parent_dir:   "create"

    working_dir:          "/scratch/xyz1234/my_project"
```

Here we see the special `ym` key under which is stored the internal config of YAMLmake. To change YAMLmake internal settings from their default values we must assign new values as shown above. Any values not explicitly set using a `config` pipeline item will remain at their default values. See the source code for a complete list of internal YAMLmake configuration keys.

### module

Much like `include` a `module` pipeline item refers to another YAML file, which gets executed when the pipeline reaches the `module` item. But rather than inserting the YAML items at the location of the `module`, the file gets executed as a subordinate process (but still serially, not as a fork into a parallel process). This means that any configuration changes made within the module are forgotten once complete, and the remainder of the top-level pipeline continues to execute with the exact same configuration settings as if the module had never been loaded. The only result will therefore be any new output files created by actions within the module.

### action

An `action` executes commands that process your data files. The required fields of an action are `name`, `input`, `output` and `shell`. The name field is required and must not contain any spaces or special characters as it forms the basis for various log filenames and GridEngine job names.

The input field specifies one or more files that must be present in order for the action to be runnable. The output field specifies all the output files the action is guaranteed to create. If any have not been created by the time the action completes YAMLmake therefore assumes the action has failed and acts accordingly. The action will also not run if it looks like the output files have already been created and no input file has changed since then. The shell field contains the bash command that actually carries out the action using placeholders to insert the paths of the required input and output files.

The input and/or output fields can contain one or more simple named subfields, or be an arbitrary hierachy of fields containing lists and dictionaries with any number of named file paths. As a simple example, supposed there are two input files and two output files:

```
- action:
    name: "zip-test"
    input:
      first_input: "file1"
      second_input: "file2"
    output:
      results: "all.zip"
      aux_file: "stdout"
    shell: |
      echo Running {%name}
      zip {%results} {%first_input} {%second_input} > {%aux_file}
```

Above we see that, although the file path field names are nested under the input and output fields, for naming purposes they are considered to be directly under the action name space, and we refer to them directly within the shell using `{%first_input}` etc and not `{%input/first_input}` (see below for how to access nested variables within the configuration hierachy). Effectively the "input" and "output" field labels disappear from the name space and their subfields get promoted one level. This applies only to the input and output fields, and saves a lot of typing in the shell command, but also means we must avoid name collisions with other config variables due to this promotion.

To operate on sets of multiple files the input and output items can also contain any combination of four other special types of placeholder which expand the file paths in various ways, as shown in this table:

| placeholder type       | expansion type | number of jobs |
| ---                    | ---     | --- |
| `{%my_var}`            | refer to the variable `my_var`, which may be a config variable or a variable defined in the action, or an input/outfile file | no effect |
| `{*my_capture}`        | treat the overall string as a glob, insert the `*` wildcard here, writing the wildcard text into the variable `my_capture` | Spawns a separate instance of the shell command for each file matching the glob |
| `{=my_existing_list}`  | expand using each element of `my_existing_list` list variable | Spawns a separate instance of the shell command for each list element |
| `{+my_captures}`       | treat the overall string as a glob, insert the `*` wildcard here, writing all wildcard values as a list into the variable `my_captures`  | Spawn one instance of the shell command |
| `{-my_existing_list}`  | expand using each element of `my_existing_list` list variable | Spawn one instance of the shell command |

Once all placeholders in the input and output fields have been processed and expanded into lists and/or multiple jobs they become available within the shell command to fill out its own placeholders. Note that the shell field must always be of YAML's literal block scalar type (ie `shell: |` as shown) to ensure that commands written on separate lines are not merged together into a single line. To continue a single shell command over multiple lines use the back slash character at the end of each continuing line as you would in a normal shell script file.

#### `{*...}` Spawning separate jobs for each matching filename

```
  - action:
      name: "download_datasets"
      exec: "local"
      input:
        url: "metadata/{*dataset}/url.txt"
      output:
        fastq: "data/{*dataset}/{*dataset}.fastq.gz"
      shell: |
        echo Action {%name} using execution mode {%exec}
        echo Downloading dataset {*dataset} from {%url} to {%fastq}
        wget $(cat {%url}) -O {%fastq}
```

Here the `{*dataset}` placeholder in the *url* input field means "glob (i.e. search the filesystem for) this path using a * in place of `{*dataset}` and spawn a *separate job* for each path found". These globbing placeholders only trigger a filesystem glob based on the input field(s), since the output files are assumed not to exist before the action is run. However, once a list of matching input file paths has be found the value(s) of the globbing placeholder(s) are extracted and substituted into any corresponding placeholders in the output path patterns to create the corresponding output paths for each input path.

Therefore if we start off with a set of subfolders under `metadata` named to match each of the datasets which require downloading the above action will run a separate shell command to download each url and save it in a matching folder under `data`. Any required output folders are automatically created if missing. The action will fail if the specified matching output files do not get created. Time stamps are used to check that the output files are newer than the input files.

Across the multiple jobs spawned by the action the value assigned to the `{*dataset}` placeholder is always the same across all fields within each job, therefore when each shell command runs it is guaranteed that the value will be the same across the `"metadata/{*dataset}/url.txt"` input file path, the `"data/{*dataset}/{*dataset}.fastq.gz"` output file path and the `{*dataset}` value itself. If there had been multiple input file paths per job, all containing the same placeholder, then jobs would only be created where a match was found for all the input paths containing the same value of the placeholder. 

There is also an alternative form of globbing placeholder of the form `{+dataset}` which generates a list of items within a single job. This will be explained later.

#### `{=...}` Spawning separate jobs for each item of an existing list

The special `{=...}` placeholder expands to generate a separate job for each item in an existing config list, and can appear in the input and/or output fields of an action. Lists of strings can be defined anywhere in the configuration. Additionally each action can define any temporary local configuration it needs - this does not persist after the action has ended. Therefore a short list of files can be hard coded into the action itself (although in practise it might be better kept in a separate configuration file that is included):

```
  - action:
      name: "decompress_samples"
      exec: "local"
      sample:
        - frog
        - toad
        - newt
        - caecilian
      input:
        gzip: "data/{=sample}/{=sample}.fastq.gz"
      output:
        fastq: "data/{=sample}/{=sample}.fastq"
      shell: |
        echo Decompressing data for sample {=sample}...
        gunzip --stdout {%gzip} > {%fastq}
```

The above spawns a separate job for each member of the `sample` list. If more than one such placeholder is present in a path then all combinations of both lists are used to generate separate jobs:

```
  - action:
      name: "decompress_samples"
      exec: "local"
      sample:
        - frog
        - toad
        - newt
        - caecilian
      treatment:
        - 1A
        - 1B
        - 2
        - 3
      input:
        fastq: "data/{=sample}/{=sample}.fastq"
        conf: "protocol/{=treatment}.conf"
      output:
        processed: "results/{=sample}/{=treatment}.csv"
      shell: |
        echo Analysing sample {=sample} using method {=treatment}...
        ACMEProcessTool --input-file={%fastq} \
                        --protocol-file={%conf} \
                        --results-file={%processed}
```

The above would generate 16 jobs in total from all combinations of samples and treatments.

The variable names `sample` and `treatment` are therefore used twice in the above, once as the names of lists in the local config of the action, and again as job creation placeholders referring to those existing lists. To distinguish between these two the first character of the placeholder is `%` for normal configuration variables, and `=` for list-items-to-separate-jobs placeholders. In the shell command it is therefore valid to refer to the original list as `{%sample}` and to the value assigned to the current job as `{=sample}`. The caveat here is that the `{%sample}` placeholder always still refers to the whole list of strings rather than a single string so YAMLmake will complain that it does not know how to render the list into a substitutable string form. How to do this is fully explained next, but in outline, you must either specify which list item you want, or provide a separator string to be used to join all the items together into a single string.

As with the globbing placeholder there is also an alternative form of list expansion placeholder of the form `{-sample}` which generates a list of items within a single job instead of separate jobs. This will be also explained later.

#### Accessing values stored in lists and dictionaries

YAMLmake configuration can be hierachically structured, meaning that lists (simple ordered sequences of item) and dictionaries (where items are stored under named keys) can be nested inside each other. An example is the internal YAMLmake configuration which is stored under the `ym` subkey of the global configuration. For example:

```
- config:
    ym:
      prefix: "my_project."
```

The above is a key `prefix` within a dictionary `ym` within the top-level configuration dictionary `config`. You can setup your own structured configuration if required, for example to add a new subtree of configuration use something like the following, which could be kept in a separate config file and linked to using an `include`:

```
- config:
    metadata:
      samples:
        - toad:
            n_samples: 10
            location: "rm 7"
        - frog:
            n_samples: 12
            location: "rm 9"
            note: "moved from rm 8"
        - newt:
            n_samples: 5
            location: "rm 8"
      treatments:
        - 1A
        - 1B
        - 2
        - 3 
      
```

Note: although numbers can be entered, everything is converted to a string when loaded. You may want to put quotes around every string to prevent the precision of numerical values from being altered during the loading process, as otherwise the YAML parser will initially convert numerical values into a numerical type which will then be converted into a string using Python's `str` function.

Anywhere a string can be entered as a leaf value (ie not a dictionary key / field name) in the configuration file a variable placeholder can to inserted. Where the variable required is within a nested hierachy the full path to the variable must be used, using forward-slash separators. The value of the variable with then be inserted whenever the field value is about to be used by an action. For example, to refer to the location of the newt samples you would use a place holder like:

```
{%metadata/samples/newt/location} #gives "rm 8"
```

Within a `config` pipeline item you are setting values in the global configuration namespace which persists across the pipeline steps, and placeholders referring to this begin with a % and are surrounded by curly brackets. Within an action item any variables you set there are local to that action. When the action is carried out by the pipeline the action's local config begins as a copy of the current global configuration, which then gets modified by adding the action's locally defined configuration, such that the local field values overwrite any existing values from global configuration. The input and output file paths do not yet exist when the action's local configuration is setup, so local variables cannot refer to the, but they become available in the action's shell command.

To refer to an item in a list use the 0-base integer location as the key following the list name. Negative indexes refer to items counting backwards from the end, as with Python lists:

```
{%metadata/treatments/0} # gives "1A"
{%metadata/treatments/1} # gives "1B"
{%metadata/treatments/-1} # gives "3"
{%metadata/treatments/-2} # gives "2"
```

To get a comma separated list of all the items use a comma as the key (any non-numerical character(s) used as the index to a list will generate a list of all items using that as the separator). The special index `N` yields the list's length. For dictionaries a special empty string key yields a list of all the keys, which can then be indexed like a normal list. So:

```
{%metadata/treatments/ }    # gives "1A 1B 2 3"
{%metadata/treatments/,}    # gives "1A,1B,2,3"
{%metadata/treatments/N}    # gives "4"
{%metadata/treatments/}     # gives "1A1B23"
<{%metadata/treatments/><}> # gives "<1A><1B><2><3>"
{%metadata/samples//N}      # gives "3", the number of samples
{%metadata/samples//0}      # gives "toad" 
```

The input and output sections of actions can also contain the four special types of placeholder which expand to generate either separate jobs or lists of items within a single job. The later therefore expand a single key into a list, which can be used within the shell command.

To insert environment variables from the shell used to invoke YAMLmake another type of placeholder is available, `{$...}`, which can be used anywhere. Note this gives the environment variable not of the shell command embedded in an action but of the shell used to invoke YAMLmake to run the entire pipeline. To get the environment variable value of the command the action is running use the normal `${}` form of shell variable:

```
- action:
    name: "shell-test"
    yamlmake-env1: "User {$USER} invoked YAMLmake on host {$HOSTNAME}"
    input:
      message: "data/message.txt"
    output:
      result: "test/output.txt"
    shell: |
      echo YAMLmake invocation information: {%yamlmake-env1}
      echo Shell placeholders also work in the shell section
      echo YAMLmake was invoked from directory {$PWD}
      cd data
      echo But this action is now running in ${PWD}
```

Normal shell variables of the form `${USER}`, `${HOSTNAME}` etc are simply ignored by YAMLmake as passed through to the shell command, therefore they give you the value when the command was actually executed, whereas `{$USER}` gets substituted by YAMLmake before the action's shell command is run.

#### `{+...}` Creating lists of file paths within single jobs by matching filenames

Similarly to the `{*dataset}` style placeholder there is an alternative expanding globbing placeholder of the form `{+dataset}` which can be used in the input and output sections of an action. These also expand by matching filenames in the filesystem but generate lists of file paths within the same job instead of generating separate jobs for each path. Unlike the `{*dataset}` style placeholder this means that the files paths are in a list which must be dealt with in the shell command. For example:

```
  - action:
      name: "download_datasets"
      exec: "local"
      input:
        url: "metadata/{+dataset}/url.txt"
      output:
        fastq: "data/{+dataset}/{+dataset}.fastq.gz"
      shell: |
        echo Downloading dataset list: {+dataset/,}
        echo Number of datasets: {+dataset/N}
        my_download_script.py --url_list={%url/,} \
                              --output_list={%fastq/,}
```

Assuming the metadata folder had subfolders called one, two and three which contained a file called url.txt this would generate the command:

```
    echo Downloading dataset list: one, two, three
    echo Number of datasets: 3
    my_download_script.py --url_list=metadata/one/url.txt,metadata/two/url.txt,metadata/three/url.txt \
                          --output_list=data/one/one.fastq.gz,data/two/two.fastq.gz,data/three/three.fastq.gz
```

#### `{-...}` Creating lists of file paths within single jobs from existing lists

The final special placeholder uses the form `{-sample}`, and generates lists of file paths within a single job. For example:

```
  - action:
      name: "decompress samples"
      exec: "local"
      sample:
        - frog
        - toad
        - newt
        - caecilian
      input:
        gzip: "data/{-sample}/{-sample}.fastq.gz"
      output:
        fastq: "data/{-sample}/{-sample}.fastq"
      shell: |
        echo Decompressing data for samples: "{-sample/ }"
        for finput in {%gzip/ }
        do
            foutput=${finput/fastq/fastq.gz}
            gunzip --stdout ${finput} > ${foutput}
        done
```

Here I've used only the `gzip` list of input paths in the shell command, but the expected output files must still be specified to YAMLmake so it can check they were created.

#### Running on GridEngine

To run an action using GridEngine an action must have the `exec` field set to `qsub`. This will automatically spawn an array job with as many tasks as their are jobs generated by the input and output filepath expansion step (ie the `{*...}` globbing and `{=...}` list based job spawning placeholders).

The built-in default YAMLmake config includes a set of subfields under a `qsub` field which define default resource requests for qsub jobs which can be overridden in the usual way, either by changing their values in the global configuration using `config` and/or by providing per-action value within each `action`. The qsub specific values are:

```
    qsub:
        template: 'default'        #template job script: "default" or path to your own
        log_dir:  '{%ym/log_dir}'  #log dir for qsub stdout and stderr files
        time:     '02:00:00'       #$ -l h_rt={time} maximum run time
        mem:      '4G'             #$ -l mem={mem}   maximum memory per core
        tmpfs:    '10G'            #$ -l tmpfs={tmpfs}   require tmpfs file space
        pe:       'smp'            #$ -pe {pe} {cores}   parallel execution environment
        cores:    '1'              #$ -pe {pe} {cores}   how many cores
```

As indicated in the comments above the values of time,mem,tmpfs,pe and cores are copied into the qsub job script, whereas the log_dir field can be used to set a different log directory from the global YAMLmake log directory if desired (the default, as can be seen, is to use the same). Finally, the template file used to generate the qsub job script can be overridden if you need to customise it for you local system. The default is called qsub_template.sh within mutein/pymodules. You can leave all of these `qsub` settings in the action to document the job's resource requirement even if you then decide to run it locally by setting `exec` to `local` or `parallel` (see below), and YAMLmake will just ignore them.

#### Sequential Local Execution

The default execution mode is local execution, which simply runs one job at a time without any kind of parallelism on the machine used to invoke YAMLmake. This is good for doing simple file operations that are allowed to run on the HPC login node without requiring a compute node. Select this execution mode using the `action` key `exec: "local"`. Note this execution mode does not pay any attention to the memory or CPU core settings under the `qsub` key, so it is up to you to make sure the local environment has sufficient cores, memory etc, and that you are permitted to run this type of job there.

#### Aggregated Local Execution

Sometimes you have a large number of trivial operations to perform, such as renaming a single file each time, and it would be a waste of time to wait for compute nodes to be allocated for each one, therefore you want to use local execution, but the setup time for YAMLmake to spawn each local job is making the process very slow. To speed this up it is possible to tell YAMLmake to run multiple local sequential jobs within each subshell. This has the advantage of speeding things up, at the cost of forcing multiple jobs to share the same bash session, which can cause problems in some situations. It is up to you to make sure the `shell` command part of the `action` can be run multiple time in the same shell without problems (such as ending up in the wrong working directory at the end of the first job). To enable this mode use the `local` setting of the `exec` option, add the `aggregate` option under the `ym` key, and set the value to the number of jobs you want to run in each bash subshell:

```
- action:
    name: "rename_files"
    exec: "local"
    ym:
        aggregate: "40"
    ...
```

Naturally you could just as well write your job to process multiple files within a single job from the beginning by using a list based placeholder in conjunction with a bash for loop in the shell command, but then if you need to refactor your pipeline to go back to having one file processed per job so any reason then you'd need to change the shell code, whereas with the `aggregate` method all you have to do is disable to option.

#### Parallel Local Execution

If you have direct access to a multicore machine, such as the `skinner` node on CS cluster, the you will want to run multiple jobs at the same time on the local machine, which is possible using the `parallel` execution mode. This mode takes one additional parameter nested under the `ym` key, also called `parallel`, which tells YAMLmake how many *jobs* to run at once. As with `local` mode, this execution mode ignores all the resource settings under `qsub`, so it is again up to you to make sure the machine has enough resources to run the number of jobs you requested in parallel:

```
- action:
    name:  "mark_duplicates"
    exec:  "parallel"
    conda: "gatk4"

    ym:
        parallel: "24"
    ...
```

#### Loading Simple Lists

A special type of placeholder can be used to extract data from an external, non-YAML text file. For example to load a list of values from the first column of a csv file:

```
- config:
    data_list: "{>data/my_data.csv[,C0]}"
```

Where `,C0` means to extract the 0th (0-based numbering) Column using `,` as the column separator. Use `R` to instead extract the specified row. Finally to just extract the entire file into a single variable omit the square brackets section altogether:

```
- config:
    url: "{>metadata/toad/toad_url.txt}"
```

#### Embedded Includes

A special field called `includes` can be inserted anywhere into the configuration which must be a list of one or more filepaths, which are loaded into the configuration hierachy at the location of the includes field. Unlike a top level `- include` item, files loaded using the embedded includes feature should not contain any toplevel pipeline items (which are only valid at the highest level of the pipeline), rather they must contain the exact configuration that needs to be inserted. For example:

```
- config:
    base_config:
      some_variable: "my_config_value"
      includes:
        - "auxiliary_config.yml"
```

And auxiliary_config.yml could contain:

```
another_variable: "my_other_value"
some_list:
  - "item1"
  - "item2"
```

Which would be give the same result as if the original config item has been:

```
- config:
    base_config:
      some_variable: "my_config_value"
      another_variable: "my_other_value"
      some_list:
        - "item1"
        - "item2"
```

#### Notes

Includes are relative to the config files unless they start with a "./" in which case they are relative to the working directory? Special top-level field: `run:` can be "always", "never" or "conditional". Special top-level field: `env:` will be included in the running job's bash environment variables, all values must be simple strings, no nested containers allowed here.
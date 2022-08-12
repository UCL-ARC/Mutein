# YAMLmake Workflow Tool

This markdown documents YAMLmake the new workflow tool developed during the Mutein project. This is similar to Snakemake in that a config file documents each possible action in the pipeline in terms of input and output files, the shell command(s) to be run and also allows local or HPC execution. However YAMLmake config files are written in pure YAML rather than a mixture of YAML and Snakemake exhanced Python. YAMLmake also does not create a dependency graph working backwards from a target, rather the actions are carried out sequentially in the order they appear in the file. This makes it much simpler to deal with lists of input files that are not known in advance, the action simply globs then before it runs. Whereas Snakemake requires the special checkpoint concept to handle this common situation. YAMLmake allows easy changing between local and qsub execution, whereas Snakemake requires a separate directive for each rule to make it run remotely, and these can break when a module is imported. YAMLmake tries to specify all input and output files using special placeholders which insert filenames from predefined lists or from filesystem globs executed just prior to running the action. It therefore does not need to use embedded Python functions to generate the file lists.

## Key aims

- Produce a basic workflow engine that is easy to learn and not off-putting to new data scientists learning Linux and HPC on the job as part of their broader research
- No requirement to know any scripting language beyond what bash or scripts commands you would be running manually without a workflow engine
- Support reproducibility / portability through conda
- Support local and HPC execution through GridEngine
- Support array jobs on GridEngine
- Easily rerun only the failed jobs from an array
- Easily specify sets of input and output files for each action using placeholders rather than embedded Python functions
- Flexible, hierachical, modular config files with comments to support well organised, maintainable pipelines

## Possible future features

- SLURM and Singularity / Apptainer support

## YAMLmake overview

The following outline of YAMLmake assumes some knowledge of bash.

YAMLmake runs a YAML pipeline file containing a list of *actions* from top to bottom executing each one in the order encountered, taking account of *includes* and *modules*. A configuration state is also maintained as a hierachical tree of nested lists, dictionaries and strings which is modified everytime a *config* item is encountered in the pipeline. For example:

```
  - include: "base_config.yml"
  - action:
      name: "download_reference"
      exec: "local"
      input: "metadata/GRCh38_reference_url.txt"
      output: "reference/GRCh38.fasta.gz"
      shell: |
        wget $(cat {%input}) -O {%output}
```

Each action specifies one or more input and output files. If the output files are all present and more recent than any input file then the action is not executed. Likewise if any input file is missing the action will not run. The shell field of the action contains curly bracket placeholders where YAMLmake substitutes the input and output filenames before running the command either locally or through qsub. Provided the return code is zero and all the specified output files were created the action is considered to have succeeded. If the action is considered to have failed then any output files present will be removed, either by deletion or moving to a specified recycle bin folder. They can also be flagged as stale in-place by setting the mtime to a special value.

### include

The `include` item above causes the referred to file to be inserted into the pipeline as if it had been copy-pasted in at the location of the include. The included file can contain any YAMLmake configuration, including further nested includes and modules.

### config

In the example above the file `base_config.yml` is included as the first item. This file would therefore likely contain the starting configuration options for the pipeline in the form of a `config` item, for example: 

```
- config:
    ym:
      prefix:               "my_project."
      stale_output_file:    "recycle"
      stale_output_dir:     "ignore"
      missing_parent_dir:   "create"

    working_dir:          "/scratch/xyz1234/my_project"
```

Here we see the special `ym` key underwhich is stored the internal config of YAMLmake. To change YAMLmake internal settings from their default values we must assign new values as shown above. Any values not explicitly set using a `config` pipeline item will remain at their default values. See the source code for a complete list of internal YAMLmake configuration keys.

### module

Much like `include` a `module` pipeline item refers to another YAML file, which gets executed when the pipeline reaches the `module` item. But rather than inserting the YAML items at the location of the `module`, the file gets executed as a subordinate process (but still serially, not as a fork into a parallel process). This means that any configuration changes made within the module are forgotten once complete, and the remainder of the top-level pipeline continues to execute with the exact same configuration settings as if the module had never been loaded. The only result will therefore be any new output files created by actions within the module.

### action

The `input` and `output` keys of an action are special. In their simplest form they can contain a string specifying a single input and output file respectively, referred to simply as `{%input}` and `{%output}` in the shell command. However they can also be expanded into dictionaries containing any number of keys, each of which is then available within the shell command. For example:

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
      zip {%results} {%first_input} {%second_input} > {%aux_file}
```

To operate on sets of multiple files the input and output items can also contain any combination of four special types placeholders which expand the filenames in various ways. These four placeholder types are of the form `{*placeholder}`, `{=placeholder}`, `{+placeholder}` and `{-placeholder}` which respectively operate to expand the input and output lists by globbing paths that spawn separate jobs (`{*...}`), using existing variable lists to spawn separate jobs (`{=...}`), globbing paths that create path lists within single jobs (`{+...}`) and using existing variable lists to create path lists within single jobs (`{-...}`). This will be explained in more detail below.

#### Spawning separate jobs by matching filenames

```
  - action:
      name: "download_datasets"
      exec: "local"
      input:
        url: "metadata/{*dataset}/url.txt"
      output:
        fastq: "data/{*dataset}/{*dataset}.fastq.gz"
      shell: |
        echo Downloading dataset {*dataset}
        wget $(cat {%url}) -O {%fastq}
```

Here the `{*dataset}` placeholder in the *url* input key means "glob (i.e. search the filesystem for) this path using a * in place of `{*dataset}` and spawn a *separate job* for each path. The globbing is always done using the input key, since the output files do not exist before the action is run. Once the list of matching file paths has be created the value of the placeholder is extracted in a variable and substituted into all output path patterns to create the corresponding output paths for each input path.

Therefore if we start off with a set of subfolders under `metadata` named to match each of the datasets that require downloading the above action will run a separate shell command to download each url and save it in a matching folder under `data`. Any required output folders are automatically created if missing. The action will fail if the specified matching output files do not get created. Time stamps are used to check that the output files are newer than the input files.

Across the multiple jobs spawned by the action the value assigned to the `{*dataset}` placeholder is always the same across all keys within each job, therefore when each shell command runs it is guaranteed that the value will be the same across the `"metadata/{*dataset}/url.txt"` input file path, the `"data/{*dataset}/{*dataset}.fastq.gz"` output file path and the `{*dataset}` value itself.

There is also an alternative form of globbing placeholder of the form `{+dataset}` which generates a list of items within a single job. This will be explained later.

#### Spawning separate jobs from an existing list

Lists of strings can be defined anywhere in the configuration. Additionally each action can define any temporary local configuration it needs - this does not persist after the action has ended. Therefore a short list of files can be hard coded into the action itself:

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
        gzip: "data/{=sample}/{=sample}.fastq.gz"
      output:
        fastq: "data/{=sample}/{=sample}.fastq"
      shell: |
        echo Decompressing data for sample {=sample}...
        gunzip --stdout {%gzip} > {%fastq}
```

The above spawns a separate job for each member of the `sample` list. If more than one such placeholder is present in a path then all combinations of both lists are generated as separate jobs:

```
  - action:
      name: "decompress samples"
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

The variable names `sample` and `treatment` are therefore used twice in the above, once as the names of lists in the local config of the action, and again as job creation placeholders referring to the existing lists. To distinguish between these two the first character of the placeholder is `%` for normal variables, and `=` for list-to-separate-job placeholders. In the shell command it is therefore valid to refer to the original list as `{%sample}` and to the value assigned to the current job as `{=sample}`. The caveat here is that the `{%sample}` placeholder refers to a list of strings rather than a single string so YAMLmake will complain that it does not know how to render the list into a substitutable string form. How to do this is fully explained next, but in outline, you must either specify which list item you want, or provide a separator string to be used to join all the items together into a single string.

As with the globbing placeholder there is also an alternative form of list expansion placeholder of the form `{-sample}` which generates a list of items within a single job instead of separate jobs. This will be also explained later.

#### Accessing values stored in lists and dictionaries

YAMLmake configuration can be hierachically structured, meaning that lists (simple ordered sequences of item) and dictionaries (where items are stored under named keys) can be nested inside each other. An example is the internal YAMLmake configuration which is stored under the `ym` subkey of the global configuration. For example:

```
- config:
    ym:
      prefix: "my_project."
```

The above is a key `prefix` within a dictionary `ym` within the top-level configuration dictionary `config`. Youc an setup your own structured configuration if required, for example to add a new subtree of configuration use something like the following, which could be kept in a separate config file and linked to using an `include`:

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

Note: although numbers can be entered, everything is converted to a string when loaded. You may want to put quotes around everything to prevent the precision of numerical values from being altered during the loading process.

Anywhere a string can be entered into the configuration file a variable placeholder can to inserted. Where the variable required is within a nest hierachy the full path to the variable must be used, using forward-slash separators. For example, to refer to the location of the newt samples you would use a place holder like:

```
{%metadata/samples/newt/location} #give "rm 8"
```

Note: you never have to put `config` or `action` as the first part of the key. To refer to an item in a list use the 0-base integer location as the key following the list name. Negative indexes refer to items counting backwards from the end, as with Python lists:

```
{%metadata/treatments/0} # gives "1A"
{%metadata/treatments/1} # gives "1B"
{%metadata/treatments/-1} # gives "3"
{%metadata/treatments/-2} # gives "2"
```

To get a comma separated list of all the items use a comma as the key (any character(s) used as the index will generate a list of all item using that as the separator). The special index `N` yields the list's length. For dictionaries a special empty string key yields a list of all the keys, which can then be indexed like a normal list. So:

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

#### Creating lists of file paths within single jobs by matching filenames
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
        echo Downloading dataset list: {+dataset/, }
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

#### Creating lists of file paths within single jobs from existing lists
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
        echo Decompressing data for samples: "{-sample/" "}"
        for finput in {%gzip/ }
        do
            foutput=${finput/fastq/fastq.gz}
            gunzip --stdout ${finput} > ${foutput}
        done
```

Here I've used only the `gzip` list of input paths in the shell command, but the expected output files must still be specified to YAMLmake so it can check they were created.
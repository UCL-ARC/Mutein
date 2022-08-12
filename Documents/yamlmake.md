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

To operate on sets of multiple files the input and output items can contain named subitems, each of which can contain special placeholders which expand the filenames by globbing against the filesystem just before the action is run.

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

#### Spawning separate jobs by matching filenames

Here the `{*dataset}` placeholder in the *url* input key means "glob (i.e. search the filesystem for) this path using a * in place of `{*dataset}` and spawn a *separate job* for each path. The globbing is always done using the input key, since the output files do not exist before the action is run. Therefore if we start off with a set of subfolders under `metadata` named to match each of the datasets that require downloading the above action will run a separate shell command to download each url and save it in a matching folder under `data`. Any required output folders are automatically created if missing. The action will fail if the specified matching output files do not get created. Time stamps are used the check that the output files are newer than the input files.

Across the multiple jobs spawned by the action the value assigned to the `{*dataset}` placeholder is always the same across all keys within each job, therefore when each shell command runs it is guaranteed that the value will be the same across the `"metadata/{*dataset}/url.txt"` input file path, the `"data/{*dataset}/{*dataset}.fastq.gz"` output file path and the `{*dataset}` value itself.

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
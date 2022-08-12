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

To operate on sets of multiple files the input and output items can contain multiple named subitems, each of which can contain special placeholders which expand the filenames using predefined lists from the configuration object or by globbing against the filesystem just before the action is run.

```
  - action:
      name: "download_datasets"
      exec: "local"
      input:
        url: "metadata/{*dataset}/dataset_url.txt"
      output:
        md5: "data/{*dataset}/MD5.txt"
      shell: |
        wget $(cat {%url}) -O {*dataset}
```
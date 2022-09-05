# Overview of the variant calling part of the pipeline

This document explains what each of the subfolders of the main Pipeline folder contains. The subfolders currently contain multiple, incomplete versions of the same pipeline steps. The immediate aim is to allow Jeremy to find the parts he needs to begin developing on the Snakemake version of the pipeline while Rob completes the YAMLmake tool and accesses whether to port the Snakemake version over to using YAMLmake.

The pipeline is split into two main parts: the protein structure part (coded by Rachel) and the variant calling part (coded by Rob thus far, which Jeremy will now also be working on).

The foldx and geneprot folders contain the protein structure calculation code.

The basic steps of the variant calling part of the pipelines are:

- download reference genome (run on the login node as not resource intensive)
- download datasets (keogh and yokoyama only so far, run on the login node as not resource intensive)
- fastqc quality assessment of raw data (run on HPC compute nodes)
- preprocessing of the raw reads (adaptor trim, quality trim, recalibrate base quality scores etc, run on HPC compute nodes)
- mapping reads to reference (run on HPC compute nodes)
- calling germline variant (run on HPC compute nodes)
- calling somatic variants (run on HPC compute nodes, not yet implemented)

Some of the files and folders relating to the pipeline versions no longer under development have now been moved into a new subfolder under Pipelines called "old". This may have broken some path references in the old pipeline scripts.

Variant calling pipeline versions in chronological order of creation:

## Plain bash
This version used a plain bash "generator" scripts to glob the list of data files to be processed and output corresponding lists of commands to text files, which were then executed using qsub bash scripts. The drawback of this method was that the qsub scripts read their main data processing commands from the text file, and therefore the script itself was just piping output from this file into a bash subshell, giving no indication of what command it was actually running. Alternative upgraded ("Varcall") versions of some of the steps of the pipeline (fastqc to read mapping) were then made along side them in the same folder. The current final step in the pipeline (calling germline mutations) was not upgraded. The scripts to download the genome reference and the datasets don't run on the compute nodes and are just simple bash scripts, so these also didn't change when moving to the Varcall version of the pipeline.

Folders/files:
- old/variant_calling/*
- scripts/haplotype_caller_generator.sh
- premapping/
  - get_dataset*.sh
  - generate_testdata_from_keogh2018.sh
  - get_reference.sh
- mapping/
  - bwa_index.sh
  - bwa_generator.sh
  - bwa_runner_step1.sh
  - bwa_runner_step2.sh

## Varcall
This version used python "generator" functions which output lists of input filenames to a task manifest file in csv format. The generator scripts are python and do an "import varcall as vc" at the top. Where multiple files are processed these are selecting using a function vc.glob_dirs, using config from a json file from the varcall/config folder. Jobs are then run on the HPC using the command line tool vc_submit_arrayjob, which generates a qsub script from templates in the varcall_templates folder which contain a line of bash like "$(vc_extract_params <parameter_name(s)>)"  which extracts the correct filenames from the task manifest file into bash variables for the rest of the script to use.

Folders/files:
- premapping/fastqc_generator.py
- mapping/bwa_mem_generator.py
- mapping/samtools_index_generator.py
- old/varcall_config/
  - fastqc.json
  - bwa_index.json
  - bwa_mem.json
  - dataset_keogh.json
  - etc
- old/varcall_templates
  - fastqc.qsub
  - bwa_index.qsub
  - bwa_mem.qsub

## Snakemake
The snakemake version of the pipeline only currently covers setting up conda environments, downloading and indexing the reference genome, and downloading two of the datasets (keogh and yokoyama). The conda environments used by snakemake are defined as plain lists of packages in the conda_envs folder (this is a work around as I couldn't get the official snakemake env files to work). The main snakemake config files are in the config folder. The main pipeline file is config/Snakefile.

- conda_envs
- config/
  - Snakefile

## YAMLmake
I will fully document it elsewhere, but in brief: YAMLmake uses a YAML config file format similar to Snakemake's except that no embedded python code is allowed, i.e. it is pure YAML. The file format is a YAML list at the top level, with each list item being a dictionary with a single key defining what type of item it is: "config" modifies the config, "action" defines a data processing step similar to a Snakemake rule. YAMLmake reads the file one item at a time and processes things in the order encountered (ie no complex dependency graph working backwards from a desired target). Therefore "config" items that modify the configuration affect everything after them in the file, and "actions" are processed immediately.

Configuration is arbitrary nested dicts and lists (like JSON/YAML) with strings as the leaf nodes, and addressed using path-like keys eg you could have "datasets/keogh/number_of_samples". Where the container is a list you use the python style list integer offset as the key eg "datasets/keogh/samples/-1" to access the last item in a list called samples in a dict called keogh in a dict called datasets etc. Non integer keys to lists are used as separators to the str.join() function, except "N" which returns the list length. Where the container is a dict an empty string key will return the list of keys

Actions can contain their own arbitrary config that only affects that action (not subsequent actions) and have three main special keys: input, output and shell, which are very similar to Snakemake: input contains lists of files that must all be present for the action to run and output contains the corresponding lists of output files that should be generated for the job to be considered to have succeeded.

Anywhere in the config/action you can use special placeholders similar to python fstrings that refer to other config variables or, within action input/output/shell variables, cause YAMLmake to generate lists of filenames using either predefined lists of names or using globbing, to assist in creating the list of input and output files that are then used in the shell command. These filename lists either generate separate qsub array jobs or else remain as filename lists within single jobs, or a combination of both.

Summary of YAMLmake placeholders:

- {$ENVIRONMENT_VARIABLE} "run on {$HOSTNAME} as {$USER}"
- {%YAMLMAKE_VARIABLE}
  - "the sample name is {%sample}"
  - "the second sample name is {%sample_list/1}"
  - "the second to last sample name is {%sample_list/-2}"
  - "the comma separated sample list is {%sample_list/,}"
  - "the space separated sample list is {%sample_list/ }"
  - "the XML-like sample list is <{%sample_list/><}>"
  - "the total number of samples is {%sample_list/N}"
  - "the nested dictionary variable is {%main_dict/sub_dict/final_key}"
  - "the subdictionary keys are {%main_dict/sub_dict/}"
- {+INJOB_GLOB} a glob that generates a list of filepaths within a single job, only allowed in input/output/shell inside an action
  - some_file_list: "data/{+sample}" in the action/input term to generate the list
  - "some_command --input_list={+some_file_list/,}" to pass the comma separated list to an action/shell command
  - "for x in {+some_file_list/ } ; do echo $x >> output_list ; done" to pass a space separated list to the shell
- {*SEPARATE_JOBS_GLOB} a glob that generates one array job task per glob match, only allowed in input/output/shell inside an action
  - some_file_name: "data/{*sample}" to glob the list in action/input
  - shell: "some_command --input_file={%some_file_name} --sample={*sample}" to pass the filename and sample name to an action/shell command
- {-INJOB_LIST} uses a predefined YAMLmake list that expands an action/input path pattern into a list within a single job
  - some_file_list: "data/{-sample_list}" in the action/input term to generate the list of paths from an existing list variable called sample_list
  - "some_command --input_list={-some_file_list/,}" to pass the comma separated list to an action/shell command
- {=SEPARATE_JOBS_LIST} uses a predefined YAMLmake list that expands an action/input path pattern into a separate job per list item
  - some_file_name: "data/{=sample}" in the action/input term to generate the list of paths from an existing list variable called sample
  - "some_command --input_file={%some_file_name} --sample_name={=sample}" to pass each filename to an action/shell command


Folders forming the YAMLmake version of the pipeline:

- yamlmake
- pymodules

See YAMLmake/test_pipeline.yml as the main config file, pymodules/YAMLmake.py is the python module and scripts/YAMLmake is the "executable" command.

No pipeline steps are implemented using YAMLmake yet, but it is very similarly setup to Snakemake so that porting the Snakemake steps to it should be fairly easy.

## General utils
- premapping/generate_testdata_from_keogh2018.sh generates a small test data set from keogh

- I began work on an encrypted password vault system

- scripts folder contains anything that needs to be in the PATH for execution directly by the user on the commandline.

- capture_qacct.sh captures qacct post-run job info, useful for checking why a job failed, eg if used too much ram this should tell you, be aware that gridengine jobids loop rounds to 1 again everyso often so to avoid getting info on pther people's past jobs you must either use a unique jobname and else filter the results something extra such as by userid 
# Overview of the variant calling part of the pipeline

This document explains what each of the subfolders of the main Pipeline folder contains. The subfolders currently contain multiple, incomplete versions of the same pipeline steps. The immediate aim is to allow Jeremy to find the parts he needs to begin developing on the Snakemake version of the pipeline while Rob completes the fakemake tool and accesses whether to port the Snakemake version over to using fakemake.

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

Variant calling pipeline versions in chronological order of creation:

##Plain bash
This version used a plain bash "generator" scripts to glob the list of data files to be processed and output corresponding lists of commands to text files, which were then executed using qsub bash scripts. The drawback of this method was that the qsub scripts read their main data processing commands from the text file, and therefore the script itself was just piping output from this file into a bash subshell, giving no indication of what command it was actually running. The plain bash versions of the first steps of the pipeline (fastqc to read mapping) were then "upgraded" to the next version of the pipeline (Varcall, see below), and the only folder still present as the original plain bash version is the variant_calling folder which implements the germline variant calling step, as this has not yet been upgraded to anything. The scripts to download the genome reference and the datasets don't run on the compute nodes and are just simple bash scripts, so these also didn't change when moving to the Varcall version of the pipeline.

Folders/files:
- variant_calling/*
- scripts/haplotype_caller_generator.sh
- premapping/
  - get_dataset*.sh
  - generate_testdata_from_keogh2018.sh
  - get_reference.sh

##Varcall
This version used python "generator" functions which output lists of input filenames to a task manifest file in csv format.

Folders/files:
- premapping/fastqc_generator.py
- mapping

##Snakemake
Folders forming the Snakemake version of the pipeline:
###conda_envs
###config

##Fakemake
Folders forming the fakemake version of the pipeline:
###fakemake
###pymodules
See fakemake/test_pipeline.yml as the main config file, pymodules/fakemake.py is the python module and scripts/fakemake is the "executable" command.

No pipeline steps are implemented using fakemake yet, but it is very similarly setup to Snakemake so that porting the Snakemake steps to it should be fairly easy.

##General
premapping/generate_testdata_from_keogh2018.sh generates a small test data set from keogh

encrypted password vault

scripts folder contains anything that needs to be in the PATH for execution directly by the user on the commandline.

capture_qacct.sh captures qacct post-run job info, useful for checking why a job failed, eg if used too much ram this should tell you, be aware that gridengine jobids loop rounds to 1 again everyso often so to avoid getting info on pther people's past jobs you must either use a unique jobname and else filter the results something extra such as by userid 
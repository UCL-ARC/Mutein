# Installing and Running the Variant Calling Steps of Mutein

This document covers how to setup and use the Mutein pipeline (currently it only covers the first, variant calling part of the pipeline). You will need to get conda working in order to use the pipeline. This should then allow all the other dependencies to be installed without needing access to a root or administrator account (installing conda should also not require root).

This pipeline has been tested on Linux, and is intended to run on an HPC system due to the amount of data and computation required. snakemake is used for job control and should allow running jobs on the local machine or on the HPC's compute nodes. Currently only the qsub command is implemented for HPC support (i.e. we have only implemented support for Son of Grid Engine and not SLURM etc.)

The recommended way to install and execute is to ssh into your HPC account, clone this repository into your home folder of the HPC's filesystem and execute the snakemake commands on the HPC's login node inside a [screen](https://www.gnu.org/software/screen/) session.

The steps to install and run the pipeline follow. Install everything where you will be running the pipeline's commands from, eg if you intend to run on an HPC install to your account on the HPC not your personal laptop/desktop computer.

## install conda
If you don't already have conda installed you must install it manually as the first step. On a shared HPC environment you may already have access to some version(s) of conda through a loadable module, but this may not be the most up-to-date version. In this case you may wish to simply ignore the system-provided conda and install your own version anyway.

We recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html), follow the link to access installation instructions.

## download or clone this repository
If you have the git command available you can clone the repository, otherwise you can download the repository as a zip file and extract it. Either way put the repository in a suitable folder, such as ~/repos, that is not your intended data folder:

    mkdir -p ~/repos && cd !$

and then either:

    git clone https://github.com/UCL/Mutein.git

or:

    wget https://github.com/UCL/Mutein/archive/refs/heads/main.zip #won't work for a private repo
    unzip main.zip

#(4) create data folder somewhere (not inside the repo folder)
#(5) edit mutein bootstrap script within repo
#(6) symlink to bootstrap script from data folder
source bootstrap script

### Old Instructions Follow

    #edit software_setup/mutein_settings to symlink to it from ~/.mutein_settings.sh
    source ~/.mutein_settings.sh

    configure_software_tools.sh

    get_reference.sh
    get_dataset_keogh2018.sh
    generate_testdata_from_keogh2018.sh
    get_dataset_yokoyama2019.sh

    fastqc_generator --conf dataset_keogh.json --taskfile fastqc.tasks
    vc_submit_arrayjob --taskfile fastqc.tasks --jobname fastqc --conf fastqc.json

    bwa_index_generator --conf reference.json --conf bwa_index.json --output bwa_index_tasklist

    bwa_mem_generator --conf dataset_keogh.json --taskfile bwa_mem.tasks
    vc_submit_arrayjob --taskfile bwa_mem.tasks --jobname bwa_mem --conf reference.json --conf bwa_mem.json

    haplotype_caller_generator.sh

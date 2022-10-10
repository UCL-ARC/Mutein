# Installing and Running the Variant Calling Steps of Mutein

This document covers how to setup and use the Mutein pipeline (currently it only covers the first, variant calling part of the pipeline). You will need to get conda working in order to use the pipeline. This should then allow all the other dependencies to be installed without needing access to a root or administrator account (installing conda should also not require root).

This pipeline has been tested on Linux, and is intended to run on an HPC system due to the amount of data and computation required. snakemake is used for job control and should allow running jobs on the local machine or on the HPC's compute nodes. Currently only the qsub command is implemented for HPC support (i.e. we have only implemented support for Son of Grid Engine and not SLURM etc.)

The recommended way to install and execute is to ssh into your HPC account, clone this repository into your home folder of the HPC's filesystem and execute the snakemake commands on the HPC's login node inside a [screen](https://www.gnu.org/software/screen/) session.

The steps to install and run the pipeline follow. Install everything where you will be running the pipeline's commands from, eg if you intend to run on an HPC install to your account on the HPC not your personal laptop/desktop computer.

## Installation
### install mamba (and maybe conda)

As per the docs at https://mamba.readthedocs.io/en/latest/installation.html if you already have conda you can install mamba with this command:

```bash
conda install mamba -n base -c conda-forge
```

If you don't already have conda or mamba you can install only mamba and it will give you the conda command as well.
Installers are at https://github.com/conda-forge/miniforge#mambaforge - you will want the Linux x86_64 version.

### download or clone this repository
While the repository is private (i.e. during its development stage) you can only access the code through an authorised github account. Once your account has been given access to the repository the simplest way to do a one-off download is to log into your github account, browse to the [Mutein repository](https://github.com/UCL/Mutein), click on the green Code button, then click Download zip. Then copy this zip file onto the HPC (if required).

If you will require frequent downloading of the latest version then you need to setup ssh access to github from the place you intend to execute the pipeline, explained [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account). Then you'll be able to access it direct from the command line as if it were a public repo, explained next.

If you have the git command available you can clone the repository, otherwise you can download the repository as a zip file and extract it. Either way put the repository in a suitable folder, such as ~/repos, that is not your intended data folder.

    mkdir -p ~/repos && cd !$

and then either:

    git clone https://github.com/UCL/Mutein.git

or:

    wget https://github.com/UCL/Mutein/archive/refs/heads/main.zip
    unzip main.zip

### setup your data folder
Create a new folder to contain all the project data, making sure this is not inside the repository folder. Select an appropriate folder name and location based on your particular storage setup. For example ~/mutein_data:

    mkdir -p ~/mutein_data
    cd ~/mutein_data

When running the pipeline it is assumed you will have your main data folder we just created as your working directory.

We also create a subfolder called config:

    mkdir config

Now we need to create the boot strap script that will setup your environment ready to run the pipeline:

    cd ~/mutein_data/config
    cp ~/repos/Mutein/Pipelines/config/mutein_settings .

### customise the bootstrap script
Finally edit the mutein_settings script you just copied to set the correct paths for your repository and data folders:

    nano mutein_settings

Set MUT_DIR to the repository folder (eg ~/repo/Mutein) and MUT_DATA to the data folder (eg ~/mutein_data).

## Test Your Install
### conda
### qsub

## Run the Pipeline
If running on an HPC remotely over ssh upon first logging in enter into a screen session. This will enable you to disconnect without terminating any running snakemake commands. For example:

    screen -S mutein

will create a new screen session named mutein. Search online for a gnu screen tutorial for further information. If screen is not installed you should be able to get it through conda. Alternatives are tmux or, at a push, nohup.

### source bootstrap script
The pipeline is configured to be run from the command line from the data folder. Therefore change into the data folder and set up your environment by sourcing the bootstrap script you setup during installation:

    cd ~/mutein_data
    source ./config/mutein_settings

### Old Instructions Follow

    #edit config/mutein_settings to symlink to it from ~/.mutein_settings.sh
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

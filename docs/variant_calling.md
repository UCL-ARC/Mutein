# Variant Calling

This document walks through the variant calling pipeline using the Keogh2018 dataset as an example. It assumes you have run through the [getting started](getting_started.md) guide. To understand the details of how yamlmake, Mutein's workflow tool, works you will also want to consult the documentation [here](yamlmake.md).

It assumes the processing is taking place on the CS cluster using the `skinner` compute node which is currently available exclusively for the Mutein project. If the pipeline needs to be relocated to somewhere that requires using `qsub` to submit all compute intensive jobs to an HPC queue, then the config option `exec_mode` found in `mutein/yamlmake/base_config.yml` should be set to `qsub` rather than `parallel` to prevent these jobs running on the HPC login node.

## Overview

The main pipeline is laid out in the `mutein/yamlmake/main_pipeline.yml` yamlmake file. See comments in this file and the module and include files for details. In summary:

- `setup_conda_envs.yml` module creates all the required conda environments
- `reference.yml` downloads and prepares the reference genome
- `download/setup.yml` sets up some download folders and blocks until the user has manually downloaded some required metadata

Instructions on how to perform the manual downloads are in the `download/setup.yml` file. You can either wait until the pipeline has created the `manual_downloads` folder for you, or else preemptively create this folder yourself. The pipeline will keep quiting with an error until it finds all the required manual downloads in the expected locations.

Then the pipeline is broken into one group of modules per dataset, of which only Keogh2018 is fully implemented:

- `download/keogh2018.yml` downloads and verifies the data and creates metadata files defining individuals and samples
- `prepare/keogh2018.yml` runs fastqc, the data are already trimmed so symlinks are created making fake "trimmed" files pointing to the originals
- `mapping/keogh2018.yml` maps against the reference using bwa mem, and marks duplicates and applies base quality score recalibration using GATK
- `germline/keogh2018.yml` calls germline variants using GATK's HaplotypeCaller
- `somatic/keogh2018.yml` calls somatic variants using GATK's Mutect2

In the main workspace folder the `datasets` subfolder contains the raw data, `processed` contains the modified, processed data. Under both of these folders there is a subfolder for each dataset named after the original paper's first author and year of publication. Under that there are subfolders for data subsets (Keogh2018 has only one) and under that one subfolder for each data accession. All these folders have no normal files or other subfolders so that every folder in them represents either a data subset or accession. All other content are stored under hidden folders called `.meta` to prevent yamlmake from seeing them when it performs glob operations looking for datasets.

Final output is a single TSV file `processed/keogh2018/SRP159015/.meta/all_somatic.tsv.gz`.

The same pipeline is mostly implemented for the Martincorena2015 dataset except that somatic variant calling is missing, due to lack of ability to work out which accessions are from which samples / individuals mentioned in the paper.

Data download for all the Sanger datasets was done using a temporary FTP account setup by EGA due to pyega3 tool not working. Therefore there are two pipeline files for downloading:

- `download/sanger_ega.yml` commands to download using pyega3 tool, which did not work
- `download/sanger_ega_ftp.yml` downloads data onto myriad (files are pre-encrypted) using ncftp then moved to CS and decrypted

Note files were transferred from myriad to CS using rsync running in a `screen` session (ie not using yamlmake). Downloading was done on myriad as the CS download node does not support FTP.

### Setting up

Setup the data folder. Mutein has two large storage areas on the shared HPC file system found under `/SAN/medic/Mutein` (backed up) and `/SAN/medic/MuteinScratch` (not backed up). These are only mounted on demand so you may have to `cd` into them before they show up. Use `df -h` to see how much space is available under each of these storage areas. If you were running your jobs on the HPC gridengine queue from a login node then you would want to store your data on these in order to make it accessible from all compute nodes.

The networked file system is much slower than `skinner`'s local disks, so as we will be using `skinner` for all our compute we'll be using a local disk of `skinner` for the data folder so that it is as fast as possible. `skinner` has two large local disks mounted as `/raid6` (99 terabytes) and `/scratcha` (84 terabytes), you'll need to be setup by an administrator with a subfolder on each of these that you have permission to write to. We will now setup our main data storage folder under `/scratcha` as this is the fastest local disk, and has the best chance of keeping up with demand once we start running multiple local jobs in parallel. Below `USERNAME` should be your username on the HPC:

```
mkdir -p /scratcha/USERNAME/mutein_workspace
cd /scratcha/USERNAME/mutein_workspace
```

We should also create a folder on the shared `/SAN/medic/Mutein` network share where we can store the raw downloaded data so that it will be backed up, such that if the local disks on `skinner` fail we can at least rerun the pipeline without needing to download all the data again:

```
mkdir -p /SAN/medic/Mutein/USERNAME/mutein_workspace
```

Now we need to create a config directory and populate it with settings specific for this run of the pipeline, assuming you have already cloned the repository into `/home/USERNAME/repos/Mutein`. Here I show the use of the `nano` editor (you could also use VS Code Remote SSH extension as previously mentioned) :

```
cd /scratcha/USERNAME/mutein_workspace
mkdir -p config
cp /home/USERNAME/repos/Mutein/mutein/config/mutein_settings_cs ./config/mutein_settings
nano ./config/mutein_settings
```

Now edit the `mutein_settings` file to match the directory names you have used to setup Mutein on your account:

```
export MUT_DIR=/home/USERNAME/repos/Mutein                  #mutein repository
export MUT_DATA=/scratcha/USERNAME/mutein_workspace         #main workspace
export MUT_RAW=/SAN/medic/Mutein/USERNAME/mutein_workspace  #backed up data
```

Use `CTRL-O ENTER` to save once done editing then `CTRL-X` to quit the nano editor.

Various steps require login credentials, which are obviously not stored in the Mutein repository for security reasons. Instead you must create the required files under the config directory that you just created above, set the file view permissions so that only you can view the contents, then fill out the required usernames and passwords. Then link these files to the download step(s) so that the correct account is being used to download each dataset. Look in the download pipeline files (in the `mutein/yamlmake/download` folder) eg `yokoyama2019.yml`:

- at the top of this file the config file (`mutein/yamlmake/config/yokoyama2019.yml`) is included
- this config file defines a variable `credentials` which must point to the correct credentials file
- this variable is then used in the download command in `mutein/yamlmake/download/yokoyama2019.yml`
- in this case the command is `pyega3 -cf {%credentials} files {=subset_id}` included from `mutein/yamlmake/download/download_and_verify.yml`

pyega3 requires the username and password in the following format:

```
{"username":"a.reseacher@ucl.ac.uk","password":"1234password5678"}
```

For download operations using `ncftp` , as found in the `sanger_ega_ftp.yml` file the command is `ncftpget -f {%credentials} -RT datasets/{%dataset_id} {=subset_id}` and the credentials file must be in the following format:

```
host ftp.server.institution.ac.uk
user ftp-dropbox-name
pass 1234password5678
```

Because of problems with EGA's pyega3 download system we had to request an FTP based download, therefore the pyega3 download actions are left in as placeholders in case the system ever becomes usable in the future.

### Using screen sessions
In order to allow long running jobs to continue without interruption should you encounter network disruption or need to log out of your local computer it is best to run the pipeline inside a `screen` session ([docs](https://www.gnu.org/software/screen/)). Start a new `screen` session called "pipeline":

```
screen -S pipeline
```

Now we can bootstrap the Mutein conda environment to give access to the `yamlmake` command etc:

```
source ./config/mutein_settings
```

You should see a message like:

```
Activating the mutein_main conda environment... OK
Running the mutein script... OK
To begin try:
    yamlmake --help
Or:
    yamlmake --yaml <pipeline> --dry-run
```

For convenience we will make a symlink to the repository folder containing the yamlmake pipeline file:

```
ln -s ${MUT_DIR}/mutein/yamlmake ym
```

We can now try running the first module of the variant calling pipeline in dry run mode, therefore use the following command:

```
yamlmake --yaml ./ym/main_pipeline.yml --dry-run --module setup_conda_envs.yml
```

Which should use dry run mode, ie without actually issuing any jobs. It should also create a new folder called `yamlmake_logs` containing the messages generated by the yamlmake run.

To disconnect from the `screen` session use `CTRL-A CTRL-D`. To list available screen sessions use `screen -list`. To reconnect to the session you just created used `screen -R pipeline`.

### Running Jobs for Real

It is recommended to work through the pipeline one module at a time using the `--module` option. When not in dry run mode the pipeline will terminate at the first error, allowing you to try the fix the problem before proceeding. You can also specify individual actions or sets of actions to run using the `--run-only`, `--run-from` and `--run-to` command line options (see the [docs](yamlmake.md) for details).

yamlmake will not needlessly rerun actions that have already been completed unless something has changed in the files. It considers an action to be runnable if all the required input files are present and the output files are either missing or older than the newest input file (this is the default behaviour, although it is possible to force an action to run, or to define an action without any defined input or output files, see the [docs](yamlmake.md) for details).

yamlmake keeps no explicit record of job history recording what has been run (except the logs in the yamlmake_logs folder), everything is based on what files are present and their modification timestamps. This explains why some jobs create files using the `touch` command which don't seem to do anything - they create a record that a job was completed at a certain time. This is also why intermediate files are never simply deleted once finished with, as this will trigger the action that created them to rerun, outputting new files, incorrectly making all downstream files seem out of date. Instead, to save space, a special truncation operation is used which replaces the file with an empty placeholder without changing the name or modification time. Files output from an action that exited with an error are flagged by setting their timestamp to a special value (1980/01/01) rather than being deleted so that yamlmake recognises them as stale, but also allowing you to examine the file contents to help trouble shoot.

To run the whole of the Keogh2018 pipeline you can run through each of the separate modules listed above using the following commands one at a time:

```
yamlmake --yaml ./ym/main_pipeline.yml --module setup_conda_envs.yml
yamlmake --yaml ./ym/main_pipeline.yml --module reference.yml

#this step requires the manual metadata downloads mentioned above
yamlmake --yaml ./ym/main_pipeline.yml --module download/setup.yml 

yamlmake --yaml ./ym/main_pipeline.yml --module download/keogh2018.yml
yamlmake --yaml ./ym/main_pipeline.yml --module prepare/keogh2018.yml
yamlmake --yaml ./ym/main_pipeline.yml --module mapping/keogh2018.yml
yamlmake --yaml ./ym/main_pipeline.yml --module germline/keogh2018.yml
yamlmake --yaml ./ym/main_pipeline.yml --module somatic/keogh2018.yml
```

Each of these commands will block until completed, even if they involve submitting qsub jobs to the queue and output log messages to the screen and to log files under the `yamlmake_logs` folder. For long running steps you will likely want to disconnect the screen session and monitor progress by periodically looking at the most recent log files, eg:

```
cd yamlmake_logs
ls -lhtr | tail
tail -f <most_recent_log_file>
```

Log messages contain colour highlighting information which may not display correctly on your terminal depending on how you view the file. For example if using `less` to scroll through a long log file use the following to increase readability:

```
less -SR <log_file>
```

Log files are named `ym` followed by the timestamp of when yamlmake was first invoked, followed by a suffix. The top level log file suffix is `messages`. Per action log files end with the action's name, and qsub job log files also contain the array task number and the suffix `out` or `err` denote standard out and standard error. 

#### Notes
- output and use with muteinfold


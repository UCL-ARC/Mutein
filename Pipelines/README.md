# Working on Myriad

At the time of writing UCL's main HPC is Myriad, see (here)[https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/] for further information including how to apply for an account if you don't already have one. Thd Mutein pipeline has been designed to run on Myriad, but should be relatively easy to port of any HPC running GridEngine, (and with some minor changes, to any other HPC job submission system such as SLURM).

## Connecting to Myriad using SSH

The main pipeline is set up to be run from the command line on Myriad itself. The FoldX part of the pipeline also has a GUI but this is an optional add on. Therefore to use Mutein efficiently you will need to set up your computer to access Myriad using ssh. If you use Mac or Linux as your primary operating system you will already have a terminal program with the ssh command available. From Windows you can try using the ssh command from within powershell, or install the Ubuntu App (which relies on Window's WSL system), or use a third party app such as MobaXterm or Putty. You can also create a Linux VM and run it within Windows using VirtualBox. Once you have ssh available you will need to log into myriad using your username and password. From on site you can directly connect to myriad, from off site you must either first connect to the UCL VPN or else use a jump server (a server available from outside the firewall, currently ssh.rc.ucl.ac.uk is one option), and then connect to myriad.

To avoid having to type your password each time you want to log in you can create an ssh key using the ssh-keygen command (MobaXterm and putty have their own version of this command) and then ssh-copy-id to copy it onto myriad (and the jump server). If using Mac or Linux you can config ssh to automatically use the jump server to connect to myriad. Myriad as more than one login server for load balancing purposes, and you will be allocated one at random each time you connect. This can prevent you from reconnecting to existing "screen" sessions (see below), therefore it is best to make sure you always connect to the same login server, which can also be done using your ssh config. Here is an example ssh setup, which is saved as a file called "config" under the .ssh subfolder in your home folder:

```
    Host sshjump
        Hostname ssh.rc.ucl.ac.uk
        ForwardAgent yes

    Host myriad
        Hostname login13.myriad.rc.ucl.ac.uk
        ProxyCommand ssh sshjump -W %h:%p
```

Here an alias called myriad points to a server called login13.myriad.rc.ucl.ac.uk and always connects through an intermediate server called sshjump, which is itself an alias to ssh.rc.ucl.ac.uk. This means that "ssh myriad" will connect to the login13 server of myriad through the jump server. Provided ssh-keygen and ssh-copy-id have previously been used to setup key based login to sshjump and myriad then no password will need to be typed (although it is highly recommended to secure you private key with a passphrase).

## Copying Files to and from Myriad

If you need to copy data to and from myriad from your computer this can be done using the scp command for Mac and Linux users, or else programs such as WinSCP or MobaXterm from Windows.

## Editing Files on Myriad

For quick edits to small files you can simply log into myriad using ssh and use a command line editor such as nano or vim. For any serious editing however you will want to set up a graphical edit running on my local machine that can edit files on myriad directly without having to manually copy them back and forth. On Windows MobaXterm's built in editor can handle this efficiently for single files although it cannot be used to rename files (because it actually makes a local temporary copy of remote files and then automatically copies them back on each save). To truely edit a remote file, if using Mac or Linux you can use a program called sshfs to "mount" myriad folders on your local computer so they appear as local files to any program you run. On Ubuntu (although not on the Window's Ubuntu app) for example you can install sshfs using `sudo apt install sshfs` and then to mount your myriad home and scratch folders you could use a simple script such as:

```
  #!/bin/bash
  MYRIAD_USER=<your_myriad_username_here>
  OPTIONS="-o reconnect,ServerAliveInterval=30,ServerAliveCountMax=3"
  mkdir -p ${HOME}/myriad_home
  mkdir -p ${HOME}/myriad_scratch
  sshfs ${OPTIONS} myriad:/scratch/scratch/${MYRIAD_USER} ${HOME}/myriad_scratch
  sshfs ${OPTIONS} myriad: ${HOME}/myriad_home
```

This should make your remote folders appear under local folders called myriad_home and myriad_scratch respectively. Then, in theory, you can use any programmer's text editor to edit files directly on myriad, just by browsing inside these folders on your local file system using a normal desktop app. Some editors seem not to cope well with the increased latency of remote files. I have found the free cross platform editor VS Code (not Visual Studio itself, just the VS Code editor) to be a good choice.

If you'll be developing the git repo as well, I have found that having the repo checkout onto myriad itself and editing remotely in this way works pretty well, unless you are trying to work on a train going through a tunnel, in which case a locally checked out copy of the repo would obviously be preferable.

## Tesing and Running Jobs on Myriad

Mutein uses GridEngine to submit jobs to run on Myriad. If you are testing new code that needs to run jobs on myriad it can be very slow to test and debug if you need to wait for 10-15 minutes each time a job is submitted before being able to see if any errors were produced. For this reason, for very small operations such as creating symlinks, conda environment set up or for long running file downloads it is generally acceptable to run things directly on the login node (although you should always use `top -u <your_myriad_username>` to check how much RAM you are using. But certainly for anything requiring non-trivial RAM or compute resources you must use a compute node. Therefore rather than wait each time for a job submitted using qsub to start it is best to allocate yourself an interactive qlogin session as soon as you log into myriad, and have this running inside a "screen" session, such that if you are logged out unexpectedly for any reason you can reconnect and the qlogin session will still be running. The following shows a simple example of how to do this:

From your local machine connect to myriad using ssh (or putty/mobaxterm from windows):

```
  ssh myriad
```

Then once the myriad prompt appears, create a new screen session called eg `qlogin`. If you think a previous screen session might still be running first check using:

```
  screen -list
```

This will show you any existing screen sessions which you can connect to using the session's name:

```
  screen -R < existing_name (ignore the number) >
```

If there are no sessions already create a new one using:

```
  screen -S qlogin
```

Once inside the screen session, you can request a new qlogin session (note: if reconnecting to an existing screen session you may already be in a qlogin session that you previously created as soon as the screen session reconnects). To do this use the following, but customise the resource request to suit your intended work load:

```
qlogin -pe smp 2 -l 'mem=10G,h_rt=20:00:00' -now no
```

This requests 2 cores and 10G of memory for 20 hours. The command will wait until a compute node becomes available and then log you in to it.

Once the qlogin session has started you can temporarily disconnect the screen session, and even log out of myriad and shutdown your computer, without losing the session. To cleanly disconnect the screen session using CTRL-a CTRL-d. Then to reconnect you would use `screen -R qlogin` as shown above.

## Running Mutein on Myriad using qlogin

Assuming you have set up Mutein by downloading or cloning the repository from github and setting up a working directory with the required config subfolder you would then do the following to start running Mutein tasks (assuming your mutein data folder is called mutein_data and located under your main scratch folder):

```
cd Scratch/mutein_data
source ./config/mutein_settings
```

You should now be able to run mutein commands, for example:

```
yamlmake --yaml ${MUT_DIR}/Pipelines/yamlmake/main_pipeline.yml --dryrun | less
```

An annoyance of running inside screen is that, on my terminal at least, the mouse wheel no longer makes the terminal scroll up to show previous lines of output. Therefore you may need to pipe long outputs into less and scroll using the cursor keys instead.

## Rapid Testing using qlogin

If running yamlmake tasks from within the qlogin session if you use the "local" execution mode it will run tasks one at a time immediately on the compute node your qlogin session is running on, without having to wait for them to be schedules by GridEngine (because your qlogin session has already be allocated the resources you asked for in the qlogin command). However any qsub resource request embedded in the yamlmake configuration for the action(s) you are running will be ignored using "local" execution, therefore it is up to you to check that the action(s) have enough resources within your existing qlogin session if you run things locally in this way. It makes sense to test new code using a quick qlogin local execution and then switch the action over the qsub execution mode to access the full resources it needs for a production run. Note that all you need to do is change the execution mode from local to qsub and rerun it from your qlogin session: qsub requests can be submitted from qlogin sessions, you should not need to submit them from a myriad login node.

# Pipeline Overview

There is one Pipeline subfolder for each step in the pipeline, which should be run in
the order shown below. The convention is to issue the pipeline commands from the top level folder of your project data folder, after first having added the scripts to your path by running "source ~/.mutein_settings" or adding this line to your .bashrc. The scripts should create the required subfolders within your data folder.

#### 1) software_setup: script to setup conda environments for various required tools
#### 2) premapping: downloads the reference and datasets, and performs read quality control
#### 3) mapping: maps reads against the reference genome
#### 4) variant_calling: call variants from the mapped reads
#### 10) geneprot: This starts with I think a vcf file and ends with the protein structures
#### 11) foldx: This pipeline takes a pdb file and a list of mutations and creates xxxx




# Getting Started
This guide assumes the use of the [VS Code](https://code.visualstudio.com/) editor on Windows, but any suitable setup can be used provided you can access remote linux servers via ssh and edit remote files. [MobaXterm](https://mobaxterm.mobatek.net/) is another option. Explaining how to use ssh, VS Code or MobaXterm in detail is beyond the scope of this document, but it should provide enough pointers that a few additional internet searches should provide enough information to get stared.

In order to use an HPC to accelerate the pipelines you will need to be setup with an account on the appropriate system(s). In the UCL context this currently means [Myriad](https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/) for non-sensitive, non-GDPR data, and [CS](https://hpc.cs.ucl.ac.uk/) ([contact](https://hpc.cs.ucl.ac.uk/contact-us/)) for sensitive data. This document assumes you will be using the CS system.

Once you have used this document to get a working setup on CS you can look at the [variant_calling](variant_calling.md) document to begin using the pipeline itself.

#### Setup access to the github repository
This Mutein repository is currently private, therefore you need to have a github account with permission to access the repository in order to clone it. Then you'll need to create an ssh key [link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) on the login node of the CS cluster and upload the public key from that to your github account [link](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).


### Connecting to skinner

`skinner` is a server available to the Mutein project on the CS cluster. It must be accessed using the `ssh` command, which is usually installed by default on linux or mac, and can be setup on windows in various way, such as through powershell. To access `skinner` I have set up a `.ssh/config` file under my home folder on my windows laptop (where USERNAME is my personal username) which tells ssh how to connect to skinner:

```
Host tails
  HostName tails.cs.ucl.ac.uk
  ProxyJump USERNAME@tails.cs.ucl.ac.uk
  User USERNAME

Host comic
  HostName comic.cs.ucl.ac.uk
  ProxyJump tails
  User USERNAME

Host skinner
  User USERNAME
  HostName skinner.local
  ProxyJump comic
```

`skinner` cannot be accessed directly but rather only via the `tails` jump server and the `comic` HPC login node. I have also copied my local ssh public key onto all of these servers (using the `ssh-copy-id` command) to avoid needing to type in my password each time. Log into the `skinner` server so that the command you run will execute on that server and not on one of the shared login nodes or your local machine:

```
ssh skinner
```

#### Clone the repository into your home folder
On the CS cluster you should now have an ssh key stored in your ~/.ssh folder and the public key part should be added to your github account. On the CS cluster login node make sure you can clone the repository into your home folder:
```
cd ~
mkdir -p repos && cd repos
git clone git@github.com:UCL/Mutein.git
```

#### Editing Remote Files
Using VS Code you can directly edit the mutein files on `skinner` as required. To set this up use the `Remote-SSH` [extension](https://code.visualstudio.com/docs/remote/ssh) for VS Code. To set things up click on the green icon at the bottom left with two arrows `><`, then select `Connect to host (Remote - SSH)` then select `Configure hosts`. Now tell the extension to use the `C:\Users\USERNAME\.ssh\config` file you created above as its list of hosts. You should then find `skinner` available in the list of hosts after you select `Connect to host` the next time.

#### Install mambaforge on scratch
Here we install mamba for package management (you could use miniconda instead), and tell it to store the environment files under the mutein scratch folder. This prevents the hundreds of thousands of file it creates from cloggin up CS's backup system since our scratch space is not backed up:

```
cd /SAN/medic/MuteinScratch
mkdir mamba && cd mamba
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
chmod u+x Mambaforge-Linux-x86_64.sh
./Mambaforge-Linux-x86_64.sh
```

When asked by the installer for where to install mambaforge enter:

```
/SAN/medic/MuteinScratch/mamba/mambaforge
```

When asked if you want to run conda init say yes. Once finished exit the shell then log in again. The `conda` command should now be available. I recommend using:

```
conda config --set auto_activate_base false
```

to disable activation of the default conda environment at log in.

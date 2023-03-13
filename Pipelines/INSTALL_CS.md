# Run Pipeline on CS Cluster

#### Setup access to the github repository
Mutein repository is currently private, therefore you need to have a github account with permission to access the repository in order to clone it. Then you'll need to create an ssh key (link)[https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent] on the login node of the CS cluster and upload the public key from that to your github account (link)[https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account].

#### Clone the repository into your home folder
On the CS cluster you should now have an ssh key stored in your ~/.ssh folder and the public key part should be added to your github account. On the CS cluster login node:
```
cd ~
mkdir -p repos && cd repos
git clone git@github.com:UCL/Mutein.git
```

#### Setup you working space on scratch, link to the CS settings file
Here /SAN/medic/MuteinScratch is the non-backed up working space allocated to the project on the CS cluster's shared filesystem:
```
cd /SAN/medic/MuteinScratch
mkdir 549_mutein && cd 549_mutein
mkdir config && cd config
ln -s ~/repos/Mutein/Pipelines/config/mutein_settings_cs mutein_settings
```
The mutein_settings_cs file contains environment variables defining which directories the mutein pipeline is using. Modify these if required.

#### Install mambaforge on scratch
Here we install mamba (you could use miniconda instead), and tell it to store the environment files under the mutein scratch folder. This prevents the hundreds of thousands of file it creates from cloggin up CS's backup system since our scratch space is not backed up:
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

#### Setup the yamlmake conda environment
yamlmake is the tool that manages submitting jobs to the queue for mutein, and is part of the mutein repository. Before running yamlmake from the login node command line you must setup your environment by sourcing a script from the config folder you created. The first time it runs it will create a conda environment using mamba:
```
cd /SAN/medic/MuteinScratch/549_mutein
source ./config/mutein_settings
```

`mutein_main` should now appear in front of your command prompt indicating you have the mutein_main conda environment active and are ready to run yamlmake. Test it can run by generating the help message:

```
yamlmake --help
```

To start running the pipeline proper you pass yamlmake the main pipeline file in dry run mode. There is an environment variable ${MUT_YAML} set up containing the full path to the main pipeline file therefore you can use:

```
yamlmake --dry-run --yaml ${MUT_YAML}
```

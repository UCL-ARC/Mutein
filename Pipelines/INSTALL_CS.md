# Run Pipeline on CS Cluster

On the CS cluster login node, having setup an ssh key to github.

#### Clone the repository into your home folder
```
cd ~
mkdir repos && cd repos
git clone git@github.com:UCL/Mutein.git [or git clone https://github.com/UCL/Mutein.git]
```

#### Setup you working space on scratch, link to the CS settings file
```
cd /SAN/medic/MuteinScratch
mkdir 549_mutein && cd 549_mutein
mkdir config && cd config
ln -s ~/repos/Mutein/Pipelines/config/mutein_settings_cs mutein_settings
```

#### Install mambaforge on scratch
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

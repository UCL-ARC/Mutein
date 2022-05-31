# Pipeline Overview

There is one Pipeline subfolder for each step in the pipeline, which should be run in
the order shown below. The convention is to issue the pipeline commands from the top level folder of your project data folder, after first having added the scripts to your path by running "source ~/.mutein_settings" or adding this line to your .bashrc. The scripts should create the required subfolders within your data folder.

#### 1) software_setup: script to setup conda environments for various required tools
#### 2) premapping: downloads the reference and datasets, and performs read quality control
#### 3) mapping: maps reads against the reference genome
#### 4) variant_calling: call variants from the mapped reads
#### 10) geneprot: This starts with I think a vcf file and ends with the protein structures
#### 11) foldx: This pipeline takes a pdb file and a list of mutations and creates xxxx

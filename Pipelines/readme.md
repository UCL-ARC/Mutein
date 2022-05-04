# Pipeline Overview

One subfolder for each step in the pipeline:

#### 1) software_setup: script showing steps to setup conda environments for various required tools
#### 2) reference: script to download the reference sequence against which reads will be mapped
#### 3) datasets: one script to download each dataset that will be analysed
#### 4) quality_control: scripts to run FastQC and Trim Galore on the raw reads
#### 10) GeneProt: This starts with I think a vcf file and ends with the protein structures
#### 11) FoldX: This pipeline takes a pdb file and a list of mutations and creates xxxx

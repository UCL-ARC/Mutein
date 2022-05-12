#!/bin/bash

#copy and edit example_mutein_settings
source ~/.mutein_settings.sh

configure_software_tools.sh

get_reference.sh
get_dataset_keogh2018.sh
#generate test dataset here
fastqc_generator.sh
bwa_generator.sh
haplotype_caller_generator.sh

#!/bin/bash

#copy and edit example_mutein_settings
source ~/.mutein_settings.sh

configure_software_tools.sh

get_reference.sh
get_dataset_keogh2018.sh
generate_testdata_from_keogh2018.sh
get_dataset_yokoyama2019.sh
fastqc_generator.sh
bwa_index.sh
bwa_generator.sh
haplotype_caller_generator.sh

#!/bin/bash

#edit software_setup/mutein_settings to symlink to it from ~/.mutein_settings.sh
source ~/.mutein_settings.sh

configure_software_tools.sh

get_reference.sh
get_dataset_keogh2018.sh
generate_testdata_from_keogh2018.sh
get_dataset_yokoyama2019.sh

fastqc_generator --conf dataset.json --conf fastqc.json --output fastqc_tasklist
bwa_index_generator --conf reference.json --conf bwa_index.json --output bwa_index_tasklist
bwa_mem_generator  --conf dataset.json --conf reference.json --conf bwa_mem.json --output bwa_mem_tasklist

haplotype_caller_generator.sh

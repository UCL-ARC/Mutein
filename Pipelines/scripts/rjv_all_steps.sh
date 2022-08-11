#!/bin/bash

#edit software_setup/mutein_settings to symlink to it from ~/.mutein_settings.sh
source ~/.mutein_settings.sh

configure_software_tools.sh

get_reference.sh
get_dataset_keogh2018.sh
generate_testdata_from_keogh2018.sh
get_dataset_yokoyama2019.sh

fastqc_generator --conf dataset_keogh.json --taskfile fastqc.tasks
vc_submit_arrayjob --taskfile fastqc.tasks --jobname fastqc --conf fastqc.json

bwa_index_generator --conf reference.json --conf bwa_index.json --output bwa_index_tasklist

bwa_mem_generator --conf dataset_keogh.json --taskfile bwa_mem.tasks
vc_submit_arrayjob --taskfile bwa_mem.tasks --jobname bwa_mem --conf reference.json --conf bwa_mem.json

haplotype_caller_generator.sh

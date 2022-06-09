#!/bin/bash

#edit software_setup/mutein_settings to symlink to it from ~/.mutein_settings.sh
source ~/.mutein_settings.sh

configure_software_tools.sh

get_reference.sh
get_dataset_keogh2018.sh
generate_testdata_from_keogh2018.sh
get_dataset_yokoyama2019.sh

fastqc_generator --conf dataset_test.json --taskfile dataset_test.tasks
vc_submit_arrayjob --taskfile dataset_test.tasks --jobname fastqc --conf fastqc.json

#bwa_index_generator --conf reference.json --conf bwa_index.json --output bwa_index_tasklist
vc_submit_arrayjob --n_tasks 1 --jobname bwa_index --conf bwa_index.json --conf reference.json

bwa_mem_generator --conf dataset_test.json --taskfile bwa_mem.tasks
vc_submit_arrayjob --taskfile bwa_mem.tasks --jobname bwa_mem --conf reference.json --conf bwa_mem.json

haplotype_caller_generator.sh

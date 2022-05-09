#### bwa_generator.sh
Run this to generate the bwa_joblist file and to print the required job submission commands to the screen to be manually pasted in in order to actually submit them. The subset of datasets and accessions to be processed is controlled by the same active_datasets and active_accessions files that the other job generator scripts use.

#### bwa_index.sh
The qsub script to create the bwa index of the reference genome. This must be generated before doing the actual read mapping. Edit the embedded resource requests at the top of the script as required before running the qsub command.

#### bwa_runner.sh
The qsub script to perform the actual read mapping using bwa mem. Edit the embedded qsub resource requests prior to submitting the job array.

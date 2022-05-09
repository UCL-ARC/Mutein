#### get_reference.sh
Downloads the reference genome "analysis sets" against which the reads will be mapped

#### get_dataset_keogh2018.sh
Downloads the keogh2018 dataset. Currently requires a manual download of the accession list file prior to running!

#### fastqc_generator.sh
Generates a file called fastqc_joblist and prints the required Grid Engine command to submit an array of FastQC job to the job queue. Once the command is printed to screen it must currently be manually pasted in order to actually submit it.

Currently the selection of which datasets and accession to process is done using files within the data folders called active_datasets and active_accessions, so these can be edited before running fastqc_generator.sh to restrict to a subset of the data if required. These files are currently created manually:

datasets/active_datasets: currently contains only "keogh2018"
datasets/keogh2018/active_accessions: make this a copy of Sra_Acc_List.txt to activate all accessions

#### fastqc_runner.sh
This is the script that is passed to qsub by the command printed out by fastqc_generator. Edit the embedded resource allocation requests at the top of the script as appropriate.

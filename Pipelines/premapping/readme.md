#### get_reference.sh
Downloads the reference genome "analysis sets" against which the reads will be mapped

#### get ensembl gene list
Do this manually through a web browsers using the following steps
(the minor version number may vary as long as it matches the GRCh38 reference sequence):
https://www.ensembl.org/index.html
Biomart
CHOOSE DATABASE ==> Ensembl Genes 106
CHOOSE DATASET ==> Human Genes (GRCh38.p13)
Attributes ==> Features
GENE ==> Gene stable ID, Gene Stable ID version
     ==> Chromosome name, Gene start(bp), Gene end(bp)
     ==> Gene name, Gene synonym
Results
Export all results to ==> Compressed file (.gz), TSV
GO
Save into ensembl_gene_list subfolder

#### get_dataset_keogh2018.sh
Downloads the keogh2018 dataset. Currently requires a manual download of the accession list file prior to running!
Manually copy-pasted the neurodegeneration and cancer gene lists into metadata files from Supplementary Table 1.

#### fastqc_generator.sh
Generates a file called fastqc_joblist and prints the required Grid Engine command to submit an array of FastQC job to the job queue. Once the command is printed to screen it must currently be manually pasted in order to actually submit it.

Currently the selection of which datasets and accession to process is done using files within the data folders called active_datasets and active_accessions, so these can be edited before running fastqc_generator.sh to restrict to a subset of the data if required. These files are currently created manually:

datasets/active_datasets: currently contains only "keogh2018"
datasets/keogh2018/active_accessions: make this a copy of Sra_Acc_List.txt to activate all accessions

#### fastqc_runner.sh
This is the script that is passed to qsub by the command printed out by fastqc_generator. Edit the embedded resource allocation requests at the top of the script as appropriate.

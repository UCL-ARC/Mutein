#### get_reference.sh
Downloads the reference genome "analysis sets" against which the reads will be mapped

#### get ensembl gene list
Do this manually through a web browser using the following steps
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
Downloads the keogh2018 dataset. Currently gets the accession list using entrez-direct but the complete metadata must be downloaded manually using instructions in the script.
Manually copy-pasted the neurodegeneration and cancer gene lists into metadata files from Supplementary Table 1.

#### generate_testdata_from_keogh2018.sh
Extracts the first 1000 read pairs from one sample of keogh2018 as a test set in the test dataset folder.

#### get_dataset_yokoyama2019.sh
Downloads the yokoyama2019 dataset. Although the script can in theory download the entire dataset unaided, in practise I found that there are so many timeouts that it keeps giving up and will need restarting manually several times a day. pyega3 will automatically skip over files it's already downloaded. I have set it to download each file as a single chunk as otherwise some files never pass the MD5 check if reassembled from multiple small chunks.

#### fastqc_generator.sh
Generates a file called fastqc_joblist (by default) containing a list of job commands and prints the required Grid Engine qsub command to submit an array of FastQC jobs to the job queue. Once the command is printed to screen it must currently be manually pasted in order to actually submit it.

Currently the selection of which datasets and accession to process is done using files within the data folders called active_datasets, active_subsets and active_accessions, so these should be edited before running fastqc_generator.sh to restrict to a subset of the data if required. These files are currently created manually:

#### fastqc_runner.sh
This is the script that is passed to qsub by the command printed out by fastqc_generator. Edit the embedded resource allocation requests at the top of the script as appropriate.

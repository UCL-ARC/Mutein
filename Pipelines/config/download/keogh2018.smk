import os

dataset_id="keogh2018"
subset_id="SRP159015"
subset_dir=f"datasets/{dataset_id}/{subset_id}" #only one subset for this paper
accession_file=os.path.join(subset_dir,"accession_list")

#download list of accession names
localrules: download_accession_list
checkpoint download_accession_list:
    output:
        accession_file
    conda:
        os.environ["MUT_PREFIX"]+"enabrowsertools"
    threads:
        1
    shell:
        #grab list of sra run accessions from NCBI (US site)
        '''
        esearch -db sra -query {subset_id} \
        | efetch -format xml \
        | xtract -pattern RUN -element PRIMARY_ID \
        | awk '{{print $1}}' \
        > {accession_file}
        '''

def get_filenames(wildcards):
    'return list of the first read file for each accession listed in accession_list'
    #get the name of the accession_list file path
    #we already know it but this way we only trigger this checkpoint after its been downloaded
    out = checkpoints.download_accession_list.get(**wildcards).output[0]
    with open(out) as f:
        samples = [x.strip() for x in f.read().strip().split()]
    file_list = [f"datasets/{dataset_id}/{subset_id}/{x}/{x}_1.fastq.gz" for x in samples]
    return file_list    

#require all accessions to be downloaded
#generate list of accessions from the accession file once available
#touch the download_success file once all inputs are present
localrules: require_all_accessions
rule require_all_accessions:
    input:
        get_filenames
    output:
        touch(f"datasets/{dataset_id}/download_success")

localrules: download_accession
rule download_accession:
    input:
        accession_file
    output:
        "datasets/"+dataset_id+"/"+subset_id+"/{sample}/{sample}_1.fastq.gz"
    conda:
        os.environ["MUT_PREFIX"]+"enabrowsertools"
    threads:
        1
    resources:
        #limit to one download at a time
        download_limit=1
    shell:
        #download one accession as fastq.gz file(s) from ENA (in UK)
        '''
        enaDataGet -f fastq {output.sample}
        '''

# localrules: download_accessions
# rule download_accessions:
#     input:
#         accession_file
#     output:
#         touch("datasets/{dataset}/download_success")
#     shell:
#         #download each accession as a fastq file from ENA (in UK)
#         '''
#         for accession in $(cat "{accession_file}")
#         do
#             enaDataGet -f fastq ${{accession}}
#         done
#         '''

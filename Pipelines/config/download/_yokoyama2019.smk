import os

dataset_id="yokoyama2019"
dataset_dir=f"datasets/{dataset_id}"
subset_list=["EGAD00001004464","EGAD00001004462","EGAD00010001631"]

#require all accessions to be downloaded
#touch the download_success file once all inputs are present
localrules: require_all_accessions
rule require_all_accessions:
    input:
        expand(dataset_dir+'/{subset_id}/download_success',subset_id=subset_list)
    output:
        touch(f"{dataset_dir}/download_success")
#create gatk index
rule create_gatk_dict:
    input:
        bgz_fasta
    output:
        bgz_fasta+'.dict'
    conda:
        os.environ["MUT_PREFIX"]+"gatk4"
    threads:
        1
    params:
        time="1:00:00", mem="6G", tmpfs="10G",
    shell:
        "rm -f {bgz_fasta}.dict"
        "gatk CreateSequenceDictionary -R {bgz_fasta}"


DATA_DIR=datasets/keogh2018/SRP159015 #only one subset for this paper
ACC_FILE=accession_list

mkdir -p ${DATA_DIR}
cd "${DATA_DIR}"

#grab list of sra run accessions from NCBI (across the pond)
esearch -db sra -query SRP159015 \
| efetch -format xml \
| xtract -pattern RUN -element PRIMARY_ID \
| awk '{print $1}' \
> ${ACC_FILE}

#download each accession as a fastq file from ENA (this side of the pond)
for x in $(cat "${ACC_FILE}")
do
    enaDataGet -f fastq ${x}
done \
> download-${TIMESTAMP}.stdout \
2> download-${TIMESTAMP}.stderr

#set all samples as active (manually edit active_acc_list to work selectly)
cp ${ACC_FILE} active_acc_list

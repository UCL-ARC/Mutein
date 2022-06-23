#require the reference fasta file, and its readme and md5 file
#md5 check the reference file
localrules: check_reference
rule check_reference:
    input:
        md5all=os.path.join(config["ref"]["tmp"],config["ref"]["md5"]),
        readme=os.path.join(config["ref"]["tmp"],config["ref"]["readme"]),
        fasta=os.path.join(config["ref"]["tmp"],config["ref"]["fasta"])
    output:
        md5ref=os.path.join(config["ref"]["tmp"],"MD5.txt")
    shell:
        #extract and check just the relevant checksum
        "cat {input.md5all} | grep -e '{config[ref][fasta]}' > {output.md5ref} && "
        "cd {config[ref][tmp]} && "
        "md5sum --check MD5.txt"

#download a file from config[ref][uri] to config[ref][tmp]
localrules: wget_ref_file
rule wget_ref_file:
    output:
        os.path.join(config["ref"]["tmp"],"{filename}")
    shell:
        "wget -O {output} {config[ref][uri]}/{wildcards.filename}"

#convert gzip into an indexed bgzip format
rule create_bgzip:
    input:
        gz_fasta=os.path.join(config["ref"]["tmp"],config["ref"]["fasta"]),
        md5ref=os.path.join(config["ref"]["tmp"],"MD5.txt") #require that MD5 was valid
    output:
        bgz_fasta,
        bgz_fasta+".gzi"
    conda:
        os.environ["MUT_PREFIX"]+"bwa"
    threads:
        1
    params:
        time="1:00:00", mem="6G", tmpfs="10G",
    shell:
        "gunzip --to-stdout {input.gz_fasta} | bgzip -i -I {bgz_fasta}.gzi > {bgz_fasta}"

#create fasta file index with samtools
rule create_fai:
    input:
        bgz_fasta
    output:
        bgz_fasta+".fai"
    conda:
        os.environ["MUT_PREFIX"]+"bwa"
    threads:
        1
    params:
        time="1:00:00", mem="6G", tmpfs="10G",
    shell:
        "samtools faidx {bgz_fasta}"

#create bwa index
rule create_bwtsw:
    input:
        bgz_fasta
    output:
        bwt_indexes
    conda:
        os.environ["MUT_PREFIX"]+"bwa"
    threads:
        1
    params:
        time="2:00:00", mem="6G", tmpfs="10G",
    shell:
        "bwa index -a bwtsw {bgz_fasta}"

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
        "mutein recycle {bgz_fasta}.dict END && "
        "gatk CreateSequenceDictionary -R {bgz_fasta} -O {bgz_fasta}.dict"

import os

#define the main reference genome fasta path and all its index files
bgz_fasta = os.path.join(config["ref"]["dir"],config["ref"]["fasta"])
exts = ['gzi','fai','bwt','pac','ann','amb','sa']
all_ref_files = [ bgz_fasta ] + expand('{base}.{ext}',base=bgz_fasta,ext=exts)

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
    output: os.path.join(config["ref"]["tmp"],"{filename}")
    shell:  "wget -O {output} {config[ref][uri]}/{wildcards.filename}"

#convert gzip into bgzip format, index with samtools and bwa
rule gzip_to_bgzip:
    input:
        gz_fasta=os.path.join(config["ref"]["tmp"],config["ref"]["fasta"])
    output:
        all_ref_files
        # bgz_fasta=os.path.join(config["ref"]["dir"],config["ref"]["fasta"]),
        # fname=expand('{bgz_fasta}.{ext}',
        #     bgz_fasta=os.path.join(config["ref"]["dir"],config["ref"]["fasta"]),
        #     ext=['gzi','fai','bwt','pac','ann','amb','sa']),
    conda: "{os.environ[MUT_PREFIX]}bwa"
    threads: 1
    params:
        time="4:00:00", mem="6G", tmpfs="10G",
    shell:
        "gunzip --to-stdout {gz_fasta} | "
        "bgzip --index {bgz_fasta}.gzi > {bgz_fasta} && "
        "samtools faidx {bgz_fasta} && "
        "bwa index -a bwtsw {bgz_fasta}"



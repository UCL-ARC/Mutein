localrules: test4
rule test4:
    input:
        "test/array_test/{sample}.input",
        template="config/templates/test4.jinja2"
    output:
        "test/array_test/{sample}.output.task"
    template_engine:
        "jinja2"

localrules: aggregate_test4
rule aggregate_test4:
    input:
        "test/array_test/{sample}.output.task"

all_array_test = [ f"test/array_test/{x}.output" for x in range(1,11) ]

localrules: test5
rule test5:
    input: all_array_test

# #local test rule
# #test error handling from multiline shell fragment
# localrules: test3
# rule test3:
#     input:
#         "test/input/{sample}.fastq.gz"
#     output:
#         "test/output/{sample}.bam"
#     shell:
#         """
#         bwa mem {input} > {output}
#         """

# #local test rule
# localrules: test2
# rule test2:
#     input: "data/output/output.txt"
#     output: "data/output/test.out"
#     shell: "echo {config[ref][local_dir]} >> {output}"

# #cluster test rule using a conda env
# localrules: test
# rule test:
#     output: "data/output/output.txt"
#     input: "data/input/input.txt"
#     conda: "mutein_gatk4"
#     threads: 3
#     params: time="0:35:00", mem="3G", tmpfs="11G",
#     shell: "cat {input} > {output} && gatk >> {output}"
    
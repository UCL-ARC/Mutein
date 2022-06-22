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
    
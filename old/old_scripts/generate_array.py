
'''
this is executed from snakemake using the script directive
write command to temporary file
aggregate into a single file
spawn an array job
'''

import sys

print(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config, snakemake.params)

print(sys.argv)
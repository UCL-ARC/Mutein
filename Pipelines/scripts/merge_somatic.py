#!/usr/bin/env python3

'''
merge the separate somatic vcf files into a combined TSV type file

usage: zcat *.vcf.gz | merge_somatic.py > merged_somatic.tsv
'''

import sys

sample_names = ''

for line in sys.stdin:
    if line.startswith("##"):
        #skip vcf header
        continue

    tok = line.strip().split('\t')
    assert len(tok) == 11

    if tok[0] == "#CHROM":
        #get the two sample names from the column headers
        assert tok[8] == 'FORMAT'
        sample_names = [tok[9],tok[10]]
        continue
        
    assert sample_names != ''
    tok = sample_names + tok
    print('\t'.join(tok))
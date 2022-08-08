

from Bio import SeqIO
record_dict = SeqIO.to_dict(SeqIO.parse('chromosomes.fasta', "fasta"))

chromosome_positions = {}
with open('positions.txt') as f:
    for line in f.read().splitlines():
        if line:
            chromosome, position = line.split()
            chromosome_positions[chromosome] = int(position)


for chromosome in chromosome_positions:
    seq = record_dict[chromosome]
    position = chromosome_positions[chromosome]
    base = seq[position]
    print (chromosome, position, base)
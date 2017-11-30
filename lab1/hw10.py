#!/usr/bin/env python

def func1():
    return

seqs = []
seq = ""

read_fname = "sample_data/sample.fasta"
for line in open(read_fname):
    line = line.strip()
    if line[0] == '>':
        if seq != "":
            seqs.append([seq_id, seq])
        seq_id = line[1:]
        seq = ""
    else:
        seq += line

if seq != "":
    seqs.append([seq_id, seq])

for seq_id, seq in seqs:
    print(seq_id, seq)



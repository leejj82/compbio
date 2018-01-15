#!/usr/bin/env python
import sys
from datetime import datetime, date, time

#A=0,C=1,G=2,T=3


print >> sys.stderr, "Begin at", str(datetime.now())

min_ORF_len=200
start_codons=['ATG','GTG']
stop_codons=['TAA','TAG','TGA']

start_codons_in_num=[]
stop_codons_in_num=[]
dna_seq_in_num=[]

read_fname = "NC_002951.fna"    
for line in open(read_fname):
    line = line.strip()
    for i in range(0,len(line)):
        if line[i]=="A":
            dna_seq_in_num.append(0)
        if line[i]=="C":
            dna_seq_in_num.append(1)
        if line[i]=="G":
            dna_seq_in_num.append(2)
        if line[i]=="T":
            dna_seq_in_num.append(3)


print >> sys.stderr, "Ends at", str(datetime.now())

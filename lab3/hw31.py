#!/usr/bin/env python
import sys
from datetime import datetime, date, time

#A=0,C=1,G=2,T=3

print >> sys.stderr, "Begin at", str(datetime.now())

min_ORF_len=200
start_codons=['ATG','GTG']
stop_codons=['TAA','TAG','TGA']

def change_seq_to_num(seq):
    num=[]
    for i in xrange(len(seq)):
        if seq[i]=="A":
            num.append(0)
        if seq[i]=="C":
            num.append(1)
        if seq[i]=="G":
            num.append(2)
        if seq[i]=="T":
            num.append(3)        
    return num

start_codons_in_num=[]
stop_codons_in_num=[]
dna_seq_in_num=[]

read_fname = "NC_002951.fna"   

for line in open(read_fname):
    line = line.strip()
    dna_seq_in_num.extend(change_seq_to_num(line))

for codon in range(len(start_codons)):
    start_codons_in_num.append(change_seq_to_num(start_codons[codon]))

for codon in range(len(stop_codons)):
    stop_codons_in_num.append(change_seq_to_num(stop_codons[codon]))

def reverse_complement(list):
    reverse_complement_list=[]
    len_of_list=len(list)
    for i in range(len(list)): 
        reverse_complement_list.append(3-list[len_of_list-i-1])
    return reverse_complement_list

reversecomplement=reverse_complement(dna_seq_in_num)





print >> sys.stderr, "Ends at", str(datetime.now())

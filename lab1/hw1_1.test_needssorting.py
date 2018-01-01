#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return

seqs = []
seq = ""
numseq=[]
seq_id=[]

#a=0,c=1,g=2,t=3

read_fname = "lab01.fasta"
#read_fname = "sample_data/sample.fasta"
for line in open(read_fname):
    line = line.strip()
    if line[0] == '>':
        if seq != "":
            for i in range(0,len(seq)):
                if seq[i]=='a':
                    numseq.append(0)
                elif seq[i]=='c':
                    numseq.append(1)
                elif seq[i]=='g':
                    numseq.append(2)
                else:
                    numseq.append(3)
            seqs.append([seq_id, numseq])
        seq_id = line[1:]
        seq = ""
        numseq=[]
    else:
        seq += line

if seq != "":
    for i in range(0,len(seq)):
        if seq[i]=='a':
            numseq.append(0)
        elif seq[i]=='c':
            numseq.append(1)
        elif seq[i]=='g':
            numseq.append(2)
        else:
            numseq.append(3)
    seqs.append([seq_id, numseq])

    
lengthofseqs=len(seqs)    
lengthofnumseq=len(seqs[0][1])


def reversecomplement(arr):
    arrs = []
    i=0
    for i in range(len(arr)):
        arrs.append(3-arr[i])
    arrs=arrs[::-1]     
    return arrs


overlap_file = open("lab01.olaps.test", "w")

for k in range(1,lengthofseqs):
    compare=k
    for i in range(compare,lengthofseqs):
        for j in range(1, lengthofnumseq+1):
            if seqs[compare-1][1][lengthofnumseq-j:]==seqs[i][1][:j]:
                if j>=40:
                    print  >> overlap_file, format(compare, '03d'), format(i+1, '03d'),'F',500-j
                    break
            if seqs[compare-1][1][:j]==seqs[i][1][lengthofnumseq-j:]:
                if j>=40:
                    print >> overlap_file, format(compare, '03d'),  format(i+1, '03d'),'F',-(500-j)
                    break
            if seqs[compare-1][1][lengthofnumseq-j:]==reversecomplement(seqs[i][1])[:j]:
                if j>=40:
                    print >> overlap_file, format(compare, '03d'),  format(i+1, '03d'),'R', 500-j
                    break
            if seqs[compare-1][1][:j]==reversecomplement(seqs[i][1])[lengthofnumseq-j:]:
                if j>=40:
                    print >> overlap_file, format(compare, '03d'),  format(i+1, '03d'),'R', -(500-j)
                    break           
    print k
overlap_file.close()
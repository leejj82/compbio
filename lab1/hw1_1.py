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

tempfront=[]
temprear=[]
temp=[]
reversetemp=[]

def reversecomplement(arr):
    arrs = []
    i=0
    for i in range(len(arr)):
        arrs.append(3-arr[i])
    arrs=arrs[::-1]     
    return arrs

print >> sys.stderr, "Begin at", str(datetime.now())
overlap_file = open("lab01.olaps", "w")
for i in range(0, lengthofseqs):
    tempfront=seqs[i][1][:20]
    temprear=seqs[i][1][-20:]
    for j in range(i+1, lengthofseqs):    
        reversetemp=reversecomplement(seqs[j][1])
        for k in range(20,lengthofnumseq-20):
            temp=seqs[j][1][k:k+20]
            if temprear==temp:
                if 21+k<lengthofnumseq-20:
                    if seqs[j][1][0:k]==seqs[i][1][-20-k:-20]:
                        print >> overlap_file, seqs[i][0], seqs[j][0], "F",lengthofnumseq-k-20
            temp=reversetemp[k:k+20]        
            if temprear==temp:
                if 21+k<lengthofnumseq-20:
                    if reversetemp[0:k]==seqs[i][1][-20-k:-20]:
                        print >> overlap_file, seqs[i][0], seqs[j][0], "R",lengthofnumseq-k-20
            temp=seqs[j][1][-k-20:-k]                
            if tempfront==temp:
                if 21+k<lengthofnumseq-20:
                    if seqs[j][1][-k:]==seqs[i][1][20:20+k]:
                        print >> overlap_file, seqs[i][0], seqs[j][0], "F",k+20-lengthofnumseq
                        
            temp=reversetemp[-k-20:-k]                      
            if tempfront==temp:
                if 21+k<lengthofnumseq-20:
                    if reversetemp[-k:]==seqs[i][1][20:20+k]:
                        print >> overlap_file, seqs[i][0], seqs[j][0], "R",k+20-lengthofnumseq
                        
overlap_file.close()
print >> sys.stderr, "Ends at", str(datetime.now())
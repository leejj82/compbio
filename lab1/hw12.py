#!/usr/bin/env python

def func1():
    return

seqs = []
seq = ""
numseq=[]

#a=0,c=1,g=2,t=3

read_fname = "sample_data/sample.fasta"
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
    
def split(arr, size):
    arrs = []
    i=0
    while len(arr) > size:
         i+=1  
         pice = arr[:size]+[i]
         arrs.append(pice)
         arr   = arr[size:]
    arrs.append(arr+[i+1])
    return arrs

def my_cmp(x, y):
    assert len(x) == len(y)
    for i in range(len(x)):
        if x[i] < y[i]:
            return -1
        elif x[i] > y[i]:
            return 1
    return 0

temp=[]

for i in range(0,1):
    temp=split(seqs[i][1],10)
    temp=sorted(temp, cmp=my_cmp)


print temp
            
            






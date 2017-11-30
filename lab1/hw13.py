#!/usr/bin/env python

def func1():
    return

seqs = []
seq = ""
numseq=[]
seq_id=[]

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

temp=[]
temp1=[]
eachstrand=[]
totalarray = []

for i in range(0,10):
    temp=split(seqs[i][1],10)
    temp1=split(seqs[i][1],20)
    eachstrand.append(temp)
    eachstrand.append(temp1) 
    totalarray.append(eachstrand)
    
for i in range(0,9):
    for j in range(i+1,9):
        for l in range(0,25):
            for m in range(0,13):
                for k in range(0,len(totalarray[j][1])-9):
                    if totalarray[i][0][l][:10]==totalarray[j][1][m][k:k+10]:
                        print i,l,j,m                                





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
    
lengthofseqs=len(seqs)    

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

for i in range(0,lengthofseqs):
    temp=split(seqs[i][1],20)
    temp1=split(seqs[i][1],40)
    eachstrand.append(temp)
    eachstrand.append(temp1) 
    totalarray.append(eachstrand)
    eachstrand=[]
    
if len(seqs[0][1])%20 !=0:   
    lengthoffirstcut=len(seqs[0][1])/20+1
else:
    lengthoffirstcut=len(seqs[0][1])/20
    
if len(seqs[0][1])%40 !=0:   
    lengthofsecondcut=len(seqs[0][1])/40+1
else:
    lengthofsecondcut=len(seqs[0][1])/40

for i in range(0,lengthofseqs-1):
    for j in range(i+1,lengthofseqs-1):
        for l in range(0,lengthoffirstcut):
            for m in range(0,lengthofsecondcut):
                for k in range(0,len(totalarray[j][1][m])-19):
                    if totalarray[i][0][l][:20]==totalarray[j][1][m][k:k+20]:
                        print i+1,l,j+1,m                                

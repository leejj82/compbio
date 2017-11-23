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

temp=[]
     
for i in range(0,1):
    temp=split(seqs[i][1],10)

def my_cmp(x, y):
    assert len(x) == len(y)
    for i in range(len(x)):
        if x[i] < y[i]:
            return -1
        elif x[i] > y[i]:
            return 1
    return 0
        
def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K
    
temp=sorted(temp, key=cmp_to_key(my_cmp))

print(temp)
            
            





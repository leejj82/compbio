#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return

seqs = []
seq = ""
numseq=[]
seq_id=[]
index=0
#a=0,c=1,g=2,t=3

read_fname = "sample_data/sample.olaps"
for line in open(read_fname):
    index+=1
    line = line.strip()
    for i in range(0,len(line)):
        if line[i]==" ":
            print "yes"
        else:
            print "No"
    print index        

#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return

temp=[]
edgecollection=[]

#a=0,c=1,g=2,t=3

#read_fname = "lab01.olaps"
read_fname = "sample_data/sample.olaps"
for line in open(read_fname):
    line = line.strip()
    for i in range(0,len(line)):
        if len(temp)<2 and line[i]!=" " and  line[i+2]!=" ":
            temp.append(int(line[i:i+3])) 
        if len(temp)==2 and line[i+3]!=" ":
            temp.append(line[i+3])
        if len(temp)==3 and line[i+5]!=" ":
            temp.append(int(line[i+5:]))
    edgecollection.append(temp)
    temp=[]

def deleteoverlap(arr):
    arrs=[]
    return arrs

print deleteoverlap(edgecollection)
print edgecollection

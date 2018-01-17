#!/usr/bin/env python

import sys
import copy
from datetime import datetime, date, time

#A=0,C=1,G=2,T=3

print >> sys.stderr, "Begin at", str(datetime.now())


start_codons=['ATG','GTG']
stop_codons=['TAA','TAG','TGA']

len_start_codons=len(start_codons)
len_stop_codons=len(stop_codons)

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
    if line[0] != '>':
        dna_seq_in_num.extend(change_seq_to_num(line))
    
len_of_dna_seq_in_num=len(dna_seq_in_num)

#overlap_file = open("lab03.dnas", "w")
#for i in range(9781,10612):
#    print >> overlap_file, dna_seq_in_num[i], i+1
#overlap_file.close()




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

def compare_start_codons_with(list):
    for i in range(len_start_codons):
        if list==start_codons_in_num[i]:
            return 1
    return 0

def compare_stop_codons_with(list):
    for i in range(len_stop_codons):
        if list==stop_codons_in_num[i]:
            return 1
    return 0


def find_ORF_of_len(dna_seq):   
    all_ORF=[[],[],[]]
    #input ORF's within one strip
    for i in range(-2,len_of_dna_seq_in_num-2):
        if compare_start_codons_with(dna_seq[i:i+3])==1:
            if len(all_ORF[i%3])==0:
                all_ORF[i%3].append([i,'s'])
            elif all_ORF[i%3][-1][1]!='s':
                all_ORF[i%3].append([i,'s'])  
                        
        if compare_stop_codons_with(dna_seq[i:i+3])==1:
            if len(all_ORF[i%3])==0:
                all_ORF[i%3].append([i,'e']) 
            elif all_ORF[i%3][-1][1]!='e':
                all_ORF[i%3].append([i,'e'])  
    
    remainder=len_of_dna_seq_in_num%3
    
    for i in range(3):
        if len(all_ORF[i])!=0:
            if all_ORF[i][0][1]=='s':
                if all_ORF[(i+2)%3][-1][1]=='s':
                    all_ORF[i][0][0]=all_ORF[(i+2)%3][-1][0]-len_of_dna_seq_in_num 
                    all_ORF[(i+2)%3].pop(-1)
            if all_ORF[i][0][1]=='e':
                if all_ORF[(i+2)%3][-1][1]=='e':
                    print "Two end codons in a row without start codon. Weird"
                    all_ORF[i].pop(0)
    
    ORF_of_min_len=[[],[],[]]
    for i in range(3):
        for j in range(0,len(all_ORF[i])/2):
            if (all_ORF[i][2*j+1][0]-all_ORF[i][2*j][0]+3) >= min_ORF_len:
                ORF_of_min_len[i].append([all_ORF[i][2*j][0]+1,all_ORF[i][2*j+1][0]+3])
      
    return ORF_of_min_len
        
def maximal_ORFs(ORFlist):
    max_ORFs=[]
    for i in range(3):
        max_ORFs.extend(ORFlist[i])
    max_ORFs=sorted(max_ORFs, key=lambda ORFs: ORFs[0]) 
#    overlap_file = open("lab03.ORFs", "a")
#    for i in range(0,len(max_ORFs)):
#        print >> overlap_file, max_ORFs[i],i
#    overlap_file.close()
#    print len(max_ORFs)
    index=[]
    for j in range(len(max_ORFs)-1):
        for k in range(j+1, len(max_ORFs)):
            if max_ORFs[k][1]<max_ORFs[j][1]:
                if k not in index:
                    index.append(k)
            if max_ORFs[k][0]>max_ORFs[j][1]:
                break
    sorted(index)
#    print index
    for i in range(len(index)):
        max_ORFs.pop(index[len(index)-i-1])
#    print len(max_ORFs)
    return max_ORFs

def reverse(list):
    reverse_list=[]
    temp=[]
    for i in range(len(list)):
        reverse_list.append([len_of_dna_seq_in_num+1-list[len(list)-i-1][0],len_of_dna_seq_in_num+1-list[len(list)-i-1][1]])
    return reverse_list

def merge(ORF_of_min_len, ORF_of_min_len_reverse):
    len_P=len(ORF_of_min_len)
    len_M=len(ORF_of_min_len_reverse)
    ORFs=[]
    i,j=0,0
    while i <len_P :
        while j<len_M :
            if ORF_of_min_len_reverse[j][0] >= ORF_of_min_len[i][0]:
                ORFs.append([ORF_of_min_len[i][0],ORF_of_min_len[i][1],'+'])
                i=i+1
                if i==len_P:
                    break
            if ORF_of_min_len_reverse[j][0] < ORF_of_min_len[i][0]:
                ORFs.append([ORF_of_min_len_reverse[j][0],ORF_of_min_len_reverse[j][1],'-'])
                j=j+1
        break
    if i <len_P :
        for k in range(i, len_P):
            ORFs.append([ORF_of_min_len[k][0],ORF_of_min_len[k][1],'+'])           
    else:
        for k in range(j, len_M):
            ORFs.append([ORF_of_min_len_reverse[k][0],ORF_of_min_len_reverse[k][1],'-']) 
    return ORFs

#1
min_ORF_len=300

ORF_of_min_len=maximal_ORFs(find_ORF_of_len(dna_seq_in_num))
ORF_of_min_len_reverse=reverse(maximal_ORFs(find_ORF_of_len(reverse_complement(dna_seq_in_num))))

ORFs=merge(ORF_of_min_len, ORF_of_min_len_reverse)
        
overlap_file = open("lab03.ORF", "w")
for i in range(len(ORFs)):
    print >> overlap_file, "ORF"+format(i+1, '04d'), ORFs[i][0], ORFs[i][1], ORFs[i][2]
overlap_file.close()

#2
min_ORF_len=750

#create 61 codons
codons=[]
for i in range(4):
    for j in range(4):
        for k in range(4):
            if [i,j,k] not in stop_codons_in_num:
                codons.append([i,j,k])

#assert len(codons)==61

ORF_of_min_len=maximal_ORFs(find_ORF_of_len(dna_seq_in_num))
ORF_of_min_len_reverse=reverse(maximal_ORFs(find_ORF_of_len(reverse_complement(dna_seq_in_num))))


codoncount=[1]*61

for i in range(ORF_of_min_len):
    

overlap_file = open("lab03.ORF", "w")
for i in range(0,len(ORFs)):
    print >> overlap_file, "ORF"+format(i+1, '04d'), ORFs[i][0], ORFs[i][1], ORFs[i][2]
overlap_file.close()












print >> sys.stderr, "Ends at", str(datetime.now())

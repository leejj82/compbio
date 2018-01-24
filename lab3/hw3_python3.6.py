#!/usr/bin/env python

import sys
import copy
from datetime import datetime, date, time

#A=0,C=1,G=2,T=3

print ("Begin at", str(datetime.now()),file=sys.stderr)


start_codons=['ATG','GTG']
stop_codons=['TAA','TAG','TGA']

len_start_codons=len(start_codons)
len_stop_codons=len(stop_codons)

def change_seq_to_num(seq):
    num=[]
    for i in range(len(seq)):
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
                    print ("Two end codons in a row without start codon. Weird")
                    all_ORF[i].pop(0)
    
    ORF_of_min_len=[[],[],[]]
    for i in range(3):
        for j in range(0, int(len(all_ORF[i])/2)):
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
    print("ORF"+format(i+1, '04d'), ORFs[i][0], ORFs[i][1], ORFs[i][2], file=overlap_file)
overlap_file.close()

#2(b)
min_ORF_len=750


#create 61 codons
codons=[]#record codons
codoncount=[]#matrix of codon counts
for i in range(4):
    codoncount.append([])
    for j in range(4):
        codoncount[i].append([])
        for k in range(4):
            codoncount[i][j].append(0)
            if [i,j,k] not in stop_codons_in_num:
                codons.append([i,j,k])
                codoncount[i][j][k]=1
                

#assert len(codons)==61
#
ORF_of_min_len=maximal_ORFs(find_ORF_of_len(dna_seq_in_num))
dna_seq_in_num_rev=reverse_complement(dna_seq_in_num)
ORF_of_min_len_rev=maximal_ORFs(find_ORF_of_len(dna_seq_in_num_rev))
ORF_of_min_len_reverse=reverse(ORF_of_min_len_rev)


def codon_counting_function(ORF_of_min_len_and_rev, dna_seq_in_num_and_rev):
    for i in range(len(ORF_of_min_len_and_rev)):
        len_ORF_i=ORF_of_min_len_and_rev[i][1]-ORF_of_min_len_and_rev[i][0]+1
        for j in range(int(len_ORF_i/3)):
            assert (len_ORF_i/3)*3==len_ORF_i
            codoncount[dna_seq_in_num_and_rev[ORF_of_min_len_and_rev[i][0]-1+3*j]][dna_seq_in_num_and_rev[ORF_of_min_len_and_rev[i][0]-1+3*j+1]][dna_seq_in_num_and_rev[ORF_of_min_len_and_rev[i][0]-1+3*j+2]]+=1
            assert dna_seq_in_num_and_rev[ORF_of_min_len_and_rev[i][0]-1+len_ORF_i-3:ORF_of_min_len_and_rev[i][0]-1+len_ORF_i] in stop_codons_in_num
            
codon_counting_function(ORF_of_min_len, dna_seq_in_num)
codon_counting_function(ORF_of_min_len_rev, dna_seq_in_num_rev)

sum_of_61=0
for i in range(4):
    for j in range(4):
        for k in range(4):
            if [i,j,k] not in stop_codons_in_num:
                sum_of_61+=codoncount[i][j][k]
                
frequencies=copy.deepcopy(codoncount)

for i in range(4):
    for j in range(4):
        for k in range(4):
            if [i,j,k] not in stop_codons_in_num:
                frequencies[i][j][k]=float(frequencies[i][j][k])/sum_of_61
            else:
                frequencies[i][j][k]=0
               
overlap_file = open("lab03.expected_codons", "w")
temp="acgt"
for i in range(4):
    for j in range(4):
        for k in range(4):
            if [i,j,k] not in stop_codons_in_num:
                print( temp[i]+temp[j]+temp[k], " ", frequencies[i][j][k], file=overlap_file)
overlap_file.close()

#2(c)

tempcodoncount=[]
for i in range(4):
    tempcodoncount.append([])
    for j in range(4):
        tempcodoncount[i].append([])
        for k in range(4):
            tempcodoncount[i][j].append(0)

            
#create zero matrix
def set_to_zero(cube):
    for i in range(4):
        for j in range(4):
            for k in range(4):
                 cube[i][j][k]=0

chi_square=[]
for i in range(len(ORFs)):
    if ORFs[i][2]=='+':
        assert (ORFs[i][1]-ORFs[i][0]+1)%3==0
        codons_in_an_ORF=(ORFs[i][1]-ORFs[i][0]+1)/3-1#does not include stop codons
        for j in range(int(codons_in_an_ORF)):
            tempcodoncount[dna_seq_in_num[ORFs[i][0]+3*j-1]][dna_seq_in_num[ORFs[i][0]+3*j]][dna_seq_in_num[ORFs[i][0]+3*j+1]]+=1
                        
    else:
        codons_in_an_ORF=(ORFs[i][0]-ORFs[i][1]+1)/3-1#does not include stop codons
        for j in range(int(codons_in_an_ORF)):
            tempcodoncount[dna_seq_in_num_rev[len_of_dna_seq_in_num-ORFs[i][0]+3*j]][dna_seq_in_num_rev[len_of_dna_seq_in_num-ORFs[i][0]+3*j+1]][dna_seq_in_num_rev[len_of_dna_seq_in_num-ORFs[i][0]+3*j+2]]+=1
        
#        for k in range(4):
#            for l in range(4):
#                for m in range(4):
#                    if [k,l,m] in stop_codons_in_num:
#                        print tempcodoncount[k][l][m]
                        
    for k in range(4):
        for l in range(4):
            for m in range(4): 
                if [k,l,m] not in stop_codons_in_num:
                    tempcodoncount[k][l][m]=float(tempcodoncount[k][l][m])/codons_in_an_ORF

    chi_square.append(0)
    for k in range(4):
        for l in range(4):
            for m in range(4): 
                if [k,l,m] not in stop_codons_in_num:
                    chi_square[-1]+=(frequencies[k][l][m]-tempcodoncount[k][l][m])**2/frequencies[k][l][m]

    set_to_zero(tempcodoncount)    

for i in range(len(ORFs)):
    ORFs[i].append(chi_square[i])

for i in range(len(ORFs)):
    ORFs[i].append("ORF"+format(i+1, '04d'))
    
overlap_file = open("lab03.ORF_with_Chi_square", "w")
for i in range(len(ORFs)):
    print(ORFs[i][4], ORFs[i][0], ORFs[i][1], ORFs[i][2],ORFs[i][3], file=overlap_file)
overlap_file.close()

#we didnt consider ORFs ranging from the end to the beginning of the sequence altogether because there were none of those.

#3(b)
index=[0]*len(ORFs)
for i in range(len(ORFs)):
    if ORFs[i][2]=='+':       
        for j in range(i+1,len(ORFs)):
            if (ORFs[j][2]=='+' and ORFs[j][0]<=ORFs[i][1]) or (ORFs[j][2]=='-' and ORFs[j][1]<=ORFs[i][1]):
                if ORFs[i][3]<ORFs[j][3]:
                    index[j]+=1
                elif ORFs[i][3]>ORFs[j][3]:
                    index[i]+=1
                else:
                    print (i, 'and', j, 'have the same chi-square statistics')
#            elif ORFs[i][1]<ORFs[j][0] and ORFs[i][1]<ORFs[j][1]:# does not work because of reverse strands. Needs improvement.
#                break
    else:       
        for j in range(i+1,len(ORFs)):
            if (ORFs[j][2]=='+' and ORFs[j][0]<=ORFs[i][0]) or (ORFs[j][2]=='-' and ORFs[j][1]<=ORFs[i][0]):
                if ORFs[i][3]<ORFs[j][3]:
                    index[j]+=1
                elif ORFs[i][3]>ORFs[j][3]:
                    index[i]+=1
                else:
                    print (i, 'and', j, 'have the same chi-square statistics' )
        

    
genes=copy.deepcopy(ORFs)
 
len_index=len(index)
for i in range(len_index):
     if index[len_index-i-1]>0:
         genes.pop(len_index-i-1)
        
overlap_file = open("lab03.genes", "w")
for i in range(len(genes)):
    print(genes[i][4], genes[i][0], genes[i][1], genes[i][2],genes[i][3], file=overlap_file)
overlap_file.close()       
        
    
#3(c.1)
read_fname = "NC_002951.ptt"   

proteins=[]
for line in open(read_fname):
    line = line.split()
    proteins.append([line[0],line[1]])

proteinbase=[]
protein_stop_codon_location=[]
for protein in proteins:
    for i in range(len(protein[0])):
        if protein[0][i:i+2]=='..':
            t1=int(protein[0][0:i])
            t2=int(protein[0][i+2:])
            proteinbase.append([t1,t2,protein[1]])
            if protein[1]=='+':
                protein_stop_codon_location.append(t2)
            else:
                protein_stop_codon_location.append(t1)


gene_count=0
for i in range(len(genes)):
    stop_codon=genes[i][1]
    for j in range(len(protein_stop_codon_location)):
        if stop_codon==protein_stop_codon_location[j]:
            gene_count+=1
            break


            
##3(c.2)
sensitivity=float(gene_count)/len(proteinbase)

##3(c.3)
precision=float(gene_count)/len(genes)

overlap_file = open("lab03.counts", "w")
print("gene count", gene_count, file=overlap_file)
print("sensitivity", sensitivity,  file=overlap_file)
print( "precision", precision, file=overlap_file)
overlap_file.close()



print ("Ends at", str(datetime.now()),file=sys.stderr )

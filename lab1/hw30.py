#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return
#a=0,c=1,g=2,t=3

#numofreads=10
#readlength=250
#read_fname = "sample_data/sample.olaps"

numofreads=300 
readlength=500
read_fname = "lab01.olaps"
    

print >> sys.stderr, "Begin at", str(datetime.now())

temp=[]
edgecollection=[]
edgespernode=[] #number of edges to each 

for i in range(0, numofreads-1):
    edgespernode.append(0)


for line in open(read_fname):
    line = line.strip()
    for i in range(0,len(line)):
        if len(temp)==0 and line[i]!=" " and  line[i+2]!=" ":
            temp.append(int(line[i:i+3]))
        if len(temp)==1 and line[i+4]!=" " and  line[i+6]!=" ":
            temp.append(int(line[i+4:i+7]))
        if len(temp)==2 and line[i+8]!=" ":
            if line[i+8]=='F':
                temp.append(0)  #0 forward 1 reverse complement
            else:
                temp.append(1)
        if len(temp)==3 and line[i+10]!=" ":
            temp.append(int(line[i+10:]))
    edgecollection.append(temp)
    edgespernode[temp[0]-1]+=1
    temp=[]   

skipped=[]
temp=[]
for i in range(0,len(edgecollection)):
    if edgecollection[i][3]==0:
        skipped.append(edgecollection[i])
        for j in range(0,len(edgecollection)):
            if edgecollection[j][0]==edgecollection[i][1] or edgecollection[j][1]==edgecollection[i][1]:
                if j not in temp:
                    temp.append(j)


for j in range(0, len(temp)):
    edgespernode[edgecollection[temp[len(temp)-j-1]][0]-1]-=1  
    edgecollection.pop(temp[len(temp)-j-1])                                

nodeindex=[] #tells the starting point of each node in the array                
nodeindex.append(0)
for i in range(0,numofreads-1):
    nodeindex.append(nodeindex[i]+edgespernode[i])           
    
def findtriangles(arr):
    arrs=[]
    for i in range(0, numofreads-1): #i read number
        for j in range(nodeindex[i], nodeindex[i+1]-1): # j j nodes linked to i
            m=nodeindex[arr[j][1]-1]
            for k in range(j+1, nodeindex[i+1]): #k
                for l in range(m, nodeindex[arr[j][1]]):
                    if arr[k][1]<arr[l][1]:
                        break
                    if arr[k][1]==arr[l][1] and arr[l][0]==arr[j][1]:
                        arrs.append([j,k,l])
                        m=l+1
                        break
                    if arr[k][1]>arr[l][1]:
                        m+=1     
    return arrs

def edgestodelete(arr):
    arrs=[]
    for i in range(0, len(arr)):
        if (edgecollection[arr[i][0]][2]+edgecollection[arr[i][1]][2]+edgecollection[arr[i][2]][2])%2==0: #number of 'R' should be even
            temp=[]
            for  j in range(0,3):
                if edgecollection[arr[i][j]][3]>0 and edgecollection[arr[i][j]][2]==0:
                    temp.append([[edgecollection[arr[i][j]][0],0],[edgecollection[arr[i][j]][1],0]])
                    temp.append([[edgecollection[arr[i][j]][1],1],[edgecollection[arr[i][j]][0],1]])
                if edgecollection[arr[i][j]][3]>0 and edgecollection[arr[i][j]][2]==1:
                    temp.append([[edgecollection[arr[i][j]][0],0],[edgecollection[arr[i][j]][1],1]])
                    temp.append([[edgecollection[arr[i][j]][1],0],[edgecollection[arr[i][j]][0],1]])
                if edgecollection[arr[i][j]][3]<0 and edgecollection[arr[i][j]][2]==0:
                    temp.append([[edgecollection[arr[i][j]][1],0],[edgecollection[arr[i][j]][0],0]])
                    temp.append([[edgecollection[arr[i][j]][0],1],[edgecollection[arr[i][j]][1],1]])
                if edgecollection[arr[i][j]][3]<0 and edgecollection[arr[i][j]][2]==1:
                    temp.append([[edgecollection[arr[i][j]][1],1],[edgecollection[arr[i][j]][0],0]])
                    temp.append([[edgecollection[arr[i][j]][0],1],[edgecollection[arr[i][j]][1],0]])
            for j in range(0,6):
                for k in range(0,6):
                    for l in range(0,6):
                        if temp[j][1]==temp[k][0] and temp[l]==[temp[j][0],temp[k][1]]:
                            if arr[i][l/2] not in arrs:
                                    arrs.append(arr[i][l/2]) 
    return sorted(arrs)

triangles=findtriangles(edgecollection)

edgenumbertodelete=edgestodelete(triangles)

for i in range(0,len(edgenumbertodelete)):
    edgecollection.pop(edgenumbertodelete[len(edgenumbertodelete)-i-1])
    
alledges=[]

for i in range(0, len(edgecollection)):
    if edgecollection[i][3]>0 and edgecollection[i][2]==0:
        alledges.append([[edgecollection[i][0],0],[edgecollection[i][1],0], edgecollection[i][3]])
        alledges.append([[edgecollection[i][1],1],[edgecollection[i][0],1], edgecollection[i][3]])
    if edgecollection[i][3]>0 and edgecollection[i][2]==1:
        alledges.append([[edgecollection[i][0],0],[edgecollection[i][1],1], edgecollection[i][3]])
        alledges.append([[edgecollection[i][1],0],[edgecollection[i][0],1],edgecollection[i][3]])
    if edgecollection[i][3]<0 and edgecollection[i][2]==0:
        alledges.append([[edgecollection[i][1],0],[edgecollection[i][0],0],-edgecollection[i][3]])
        alledges.append([[edgecollection[i][0],1],[edgecollection[i][1],1],-edgecollection[i][3]])
    if edgecollection[i][3]<0 and edgecollection[i][2]==1:
        alledges.append([[edgecollection[i][1],1],[edgecollection[i][0],0],-edgecollection[i][3]])
        alledges.append([[edgecollection[i][0],1],[edgecollection[i][1],0],-edgecollection[i][3]])
        
overlap_file = open("lab01.alledges", "w")
for i in range(0,len(alledges)):
        print >> overlap_file, alledges[i]
overlap_file.close()


unitigs=[]
linkbetweenunitigs=[]

while len(alledges)!=0:
    unitig=[]
    temp=[]
    indicator=0
        
    unitig.append(alledges[0][0])
    unitig.append(alledges[0][1])
    unitig.append(alledges[0])
    alledges.pop(0)
    alledges.pop(0)
    while indicator==0: #look backward
        record=[]
        for i in range(0, len(alledges)):
            if alledges[i][1]==unitig[0]:
                record.append(i)
                for j in range(0,len(alledges)):
                    if alledges[j][0]==alledges[i][0] and j!=i:
                        record.append(j)
        if len(record)==1:
            unitig=[alledges[record[0]][0]]+[unitig[1]]+[alledges[record[0]]]+unitig[2:]
            alledges.pop(record[0])
            alledges.pop(record[0]-record[0]%2) 
            record=[]
        elif len(record)==0:
            break
        else:
            record=sorted(record)
            t=len(record)
            for i in range(0,t):
                linkbetweenunitigs.append(alledges[record[t-i-1]])
                alledges.pop(record[t-i-1])
                linkbetweenunitigs.append(alledges[record[t-i-1]-record[t-i-1]%2])
                alledges.pop(record[t-i-1]-record[t-i-1]%2)                   
            break
    while indicator==0: #look forward
        record=[]    
        for i in range(0, len(alledges)):
            if alledges[i][0]==unitig[1]:
                record.append(i)
                for j in range(0,len(alledges)):
                    if alledges[j][1]==alledges[i][1] and j!=i:
                        record.append(j)
        if len(record)==1:
            unitig.append(alledges[record[0]])
            unitig[1]=alledges[record[0]][1]
            alledges.pop(record[0])
            alledges.pop(record[0]-record[0]%2) 
            record=[]
        elif len(record)==0:
            break
        else:
            record=sorted(record)
            t=len(record)
            for i in range(0,t):
                linkbetweenunitigs.append(alledges[record[t-i-1]])
                alledges.pop(record[t-i-1])
                linkbetweenunitigs.append(alledges[record[t-i-1]-record[t-i-1]%2])
                alledges.pop(record[t-i-1]-record[t-i-1]%2)                  
            break
    unitigs.append(unitig)


overlap_file = open("lab01.unis", "w")
for i in range(0,len(unitigs)):
    sum=readlength
    for j in range(0,len(unitigs[i])-2):
        sum+=unitigs[i][j+2][2]
    print >> overlap_file, "UNI", format(i+1, '02d'), len(unitigs[i])-1, sum   
    
    print >> overlap_file, "  ", format(unitigs[i][2][0][0],'03d'), unitigs[i][2][0][1], 0  
    for j in range(1,len(unitigs[i])-2):
        print >> overlap_file, "  ", format(unitigs[i][j+2][0][0],'03d'), unitigs[i][j+2][0][1], unitigs[i][j+1][2]  
    print >> overlap_file, "  ", format(unitigs[i][len(unitigs[i])-1][1][0],'03d'), unitigs[i][len(unitigs[i])-1][1][1], unitigs[i][len(unitigs[i])-1][2]  

overlap_file.close()
print >> sys.stderr, "Ends at", str(datetime.now())


#linkbetweenunitigs
#in [i+1,k,0 or 1]
#k=0 beginning of a unitig, k=1 end of a unitig
#0=unitig left to right, 1= reverse complement

for l in range(0,2):
    for j in range(0,len(linkbetweenunitigs)):
        for i in range(0,len(unitigs)):       
            for k in range(0,2):           
                if linkbetweenunitigs[j][l]==unitigs[i][k]:
                    linkbetweenunitigs[j].append([i+1,k,0])
                if unitigs[i][k][0]==linkbetweenunitigs[j][l][0] and  unitigs[i][k][1]!=linkbetweenunitigs[j][l][1]:
                    linkbetweenunitigs[j].append([i+1,k,1])

possiblestartingpoints=[]      
for i in range(0,len(unitigs)):
    for k in range(0,2):
        index=0              
        for j in range(0, len(linkbetweenunitigs)):    
            for l in range(3,5):
                if [i+1,k]==[linkbetweenunitigs[j][l][0],linkbetweenunitigs[j][l][1]]: 
                    index=1
                    break
                    break
        if index==0:
            possiblestartingpoints.append([i+1,k])

def reversecomplement(arr):
    arrs=[]
    for i in range(0, len(arr)): 
        temp1=[arr[len(arr)-1-i][1][0],(arr[len(arr)-1-i][1][1]+1)%2]
        temp2=[arr[len(arr)-1-i][0][0],(arr[len(arr)-1-i][0][1]+1)%2]
        arrs.append([temp1,temp2,arr[len(arr)-1-i][2]])
    return arrs

import copy

def addcontig(arr1, arr2, arr3):
    arrs=copy.deepcopy(arr1)
    arrs.append([arr2[1], arr1[-1][1]+arr2[2],0,arr2[4][:2]])
    for i in range(0,len(arr3)):
        arrs.append([arr3[i][1], arrs[-1][1]+arr3[i][2],0,arr2[4][:2]])
    return arrs

#contig=[a,b,c]        
#a is a read
#b is the read position
#c=1 is if the read satisfies the condition for contig: pairs are within 2400 to 3600
       
checklist=[]
if len(possiblestartingpoints)>0:
    startingpoint=possiblestartingpoints[0]
    contig=[]
    if startingpoint[1]==0:
        contig.append([unitigs[startingpoint[0]-1][2][0],0,0,startingpoint])
        for i in range(1, len(unitigs[startingpoint[0]-1][2:])):
            contig.append([unitigs[startingpoint[0]-1][2+i][0], contig[i-1][1]+unitigs[startingpoint[0]-1][1+i][2],0,startingpoint])
        contig.append([unitigs[startingpoint[0]-1][2+i][1], contig[i][1]+unitigs[startingpoint[0]-1][2+i][2],0,startingpoint])
        
        for i in range(0, len(contig)):
            for j in range(i+1, len(contig)):
                if (contig[i][0][0]==contig[j][0][0]+1 and contig[i][0][0]%2==0) or (contig[i][0][0]==contig[j][0][0]-1 and contig[j][0][0]%2==0) :
                    if 3100>=contig[j][1]-contig[i][1]>=1900 and contig[i][0][1]!=contig[j][0][1]:
                        contig[i][2]=1
                        contig[j][2]=1
                        sum+=1
                    else: 
                        print "There is a problem in a unitig"
                    break
        for i in range(0, len(contig)):
            if contig[i][2]==0:
                checklist.append([i,contig[i][1]])
    else:
        print "Not implemented yet."
            
    index=0
    while index==0:
        linkindex=[]
        for i in range(0, len(linkbetweenunitigs)):
            if linkbetweenunitigs[i][0]==contig[-1][0]:
                linkindex.append(linkbetweenunitigs[i])
        for i in range(0, len(linkindex)):
            if linkindex[i][4][2]==0:
                temp=unitigs[linkindex[i][4][0]-1][2:]
            elif linkindex[i][4][2]==1:
                temp=reversecomplement(unitigs[linkindex[i][4][0]-1][2:])
            contigtemp=addcontig(contig,linkindex[i], temp)

        break
                
        
    
    
    
else:
    print "No possible starting points."
                

                    
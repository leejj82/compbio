#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return
#a=0,c=1,g=2,t=3

numofreads=10
read_fname = "sample_data/sample.olaps"

#numofreads=500 
#read_fname = "lab01.olaps"
    


temp=[]
edgecollection=[]
edgespernode=[] #number of edges to each 
nodeindex=[] #tells the starting point of each node in the array

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
                    if arr[k][1]==arr[l][1]:
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

edgenumbertodelete=edgestodelete(findtriangles(edgecollection))

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
    
unitigs=[]
indicator=0
    
unitigs.append(alledges[0][0])
unitigs.append(alledges[0][1])
unitigs.append(alledges[0])
alledges.pop(0)
alledges.pop(0)
record=[]
while indicator==0: #look backward
    for i in range(0, len(alledges)):
        if alledges[i][1]==unitigs[0]:
            record.append(i)
            indicator=1
            for j in range(i+1,len(alledges)):
                if alledges[j][1]==alledges[i][1]:
                    record.append(j)
                    indicator=0
            for j in range(0,len(alledges)):
                if alledges[j][0]==alledges[i][0] and j!=i:
                    record.append(j)
                    indicator=0
        if indicator==1:
            break
    if indicator==1:
        unitigs=[alledges[i][0]]+[unitigs[1]]+[alledges[i]]+unitigs[2:]
        alledges.pop(i)
        alledges.pop(i-i%2) 
        indicator=0
    else:
        break
    
while indicator==0: #look forward
    for i in range(0, len(alledges)):
        if alledges[i][0]==unitigs[1]:
            record.append(i)
            indicator=1
            for j in range(i+1,len(alledges)):
                if alledges[j][0]==alledges[i][0]:
                    record.append(j)
                    indicator=0
            for j in range(0,len(alledges)):
                if alledges[j][1]==alledges[i][1] and j!=i:
                    record.append(j)
                    indicator=0
        if indicator==1:
            break
    if indicator==1:
        unitigs.append(alledges[i])
        unitigs[1]=alledges[i][1]
        alledges.pop(i)
        alledges.pop(i-i%2) 
        indicator=0
    else:
        break


print unitigs


                


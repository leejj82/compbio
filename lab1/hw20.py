#!/usr/bin/env python
import sys
from datetime import datetime, date, time

def func1():
    return
#a=0,c=1,g=2,t=3


#numofreads=500 
#read_fname = "lab01.olaps"
    
numofreads=10
read_fname = "sample_data/sample.olaps"

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

#incomplete comparison reverse complement needed
def edgestodelete(arr):
    arrs=[]
    for i in range(0, len(arr)):
        if (edgecollection[arr[i][0]][2]+edgecollection[arr[i][1]][2]+edgecollection[arr[i][2]][2])==0: #number of 'R' should be even
            temp=[]
            for  j in range(0,3):
                if edgecollection[arr[i][j]][3]<0:
                    temp.append([edgecollection[arr[i][j]][0],edgecollection[arr[i][j]][1]])
                else:
                    temp.append([edgecollection[arr[i][j]][1],edgecollection[arr[i][j]][0]])
            for j in range(0,3):
                for k in range(0,3):
                    if j!=k and temp[j][0]==temp[k][0]:
                        for l in range(0,3):
                            if l!=j and l!=k and temp[j][1]==temp[l][1]:
                                if arr[i][j] not in arrs:
                                    arrs.append(arr[i][j])                   
        if (edgecollection[arr[i][0]][2]+edgecollection[arr[i][1]][2]+edgecollection[arr[i][2]][2])==2: #number of 'R' should be even
            print 'not yet implemented'
            
    return sorted(arrs)

edgenumbertodelete= edgestodelete(findtriangles(edgecollection))

for i in range(0,len(edgenumbertodelete)):
    edgecollection.pop(edgenumbertodelete[len(edgenumbertodelete)-i-1])
    
total=len(edgecollection)

chain=[edgecollection[0]]

for i in range(1,total):
    if chain[0][1]==edgecollection[i][0]:
        if chain[0][3]>0:
            if chain[0][2]==0:
                if edgecollection[i][2]==0:
                    if edgecollection[i][3]>0:
                        chain.append(edgecollection[i])
                        edgecollection.pop(i)
                        total-=1
                        edgecollection.append([0,0,0,0])
                    if edgecollection[i][3]<0:
                        
                if edgecollection[i][2]==1:
                    if edgecollection[i][3]>0:
                        chain.append(edgecollection[i])
                        edgecollection.pop(i)
                        total-=1
                        edgecollection.append([0,0,0,0])                    
                    if edgecollection[i][3]<0:
                        
            if chain[0][2]==1:
                if edgecollection[i][2]==0:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:
                        chain.append(edgecollection[i])
                        edgecollection.pop(i)
                        total-=1
                        edgecollection.append([0,0,0,0])                        
                if edgecollection[i][2]==1:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:
                        chain.append(edgecollection[i])
                        edgecollection.pop(i)
                        total-=1
                        edgecollection.append([0,0,0,0])
        if chain[0][3]<0:
            if chain[0][2]==0:
                if edgecollection[i][2]==0:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:
                if edgecollection[i][2]==1:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:
            
            if chain[0][2]==1:
                if edgecollection[i][2]==0:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:
                if edgecollection[i][2]==1:
                    if edgecollection[i][3]>0:
                    
                    if edgecollection[i][3]<0:


                
    if chain[0][1]==edgecollection[i][1]:
            if chain[0][2]==0:

            if chain[0][2]==1:
                









                


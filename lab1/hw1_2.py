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

overlap_file = open("lab01.triangles.test", "w")
for i in range(0,len(triangles)):
        print >> overlap_file, triangles[i] , edgecollection[triangles[i][0]],edgecollection[triangles[i][1]],edgecollection[triangles[i][2]]
overlap_file.close()

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
      
skippedlist=[]
for i in range(0, len(skipped)):
    if skipped[i][2]==0:
        skippedlist.append([[skipped[i][0],0],[skipped[i][1],0], skipped[i][3]])
        skippedlist.append([[skipped[i][1],1],[skipped[i][0],1], skipped[i][3]])
        skippedlist.append([[skipped[i][1],0],[skipped[i][0],0],-skipped[i][3]])
        skippedlist.append([[skipped[i][0],1],[skipped[i][1],1],-skipped[i][3]])
    if skipped[i][2]==1:
        skippedlist.append([[skipped[i][0],0],[skipped[i][1],1], skipped[i][3]])
        skippedlist.append([[skipped[i][1],0],[skipped[i][0],1],skipped[i][3]])
        skippedlist.append([[skipped[i][1],1],[skipped[i][0],0],-skipped[i][3]])
        skippedlist.append([[skipped[i][0],1],[skipped[i][1],0],-skipped[i][3]])
       
overlap_file = open("lab01.edges.test", "w")
for i in range(0,len(alledges)):
        print >> overlap_file, alledges[i]
overlap_file.close()

unitigs=[]

while len(alledges)!=0:
    unitig=[]
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
                alledges.pop(record[t-i-1])
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
                alledges.pop(record[t-i-1])
                alledges.pop(record[t-i-1]-record[t-i-1]%2)                   
            break
    unitigs.append(unitig)
   
index=0
for k in range(0, len(skippedlist)):#To be more precise, skipped need to be considered in all edges as well.
    for i in range(0, len(unitigs)):
        for j in range(0, len(unitigs[i])):
            if unitigs[i][j][0]==skippedlist[k][0]:
                unitigs[i][j][0]=skippedlist[k][1]
                unitigs[i].insert(j,skippedlist[k])
                index=1
            if index==1:
                break
        if index==1:
            break
    if index==1:
        break

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
                


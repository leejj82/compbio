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
contigindex=0

def contigiteration(arr1, arr4):
    linkindex=[]
    unitigtoattach=[]
    global contigprimitive
    global checklist
    for i in range(0, len(linkbetweenunitigs)):
        if linkbetweenunitigs[i][0]==arr1[-1][0]:
            linkindex.append(linkbetweenunitigs[i])
    if len(linkindex)==0:
        contigprimitive=arr1
        checklist=arr4
    for i in range(0, len(linkindex)):
        if linkindex[i][4][2]==0:
            unitigtoattach=unitigs[linkindex[i][4][0]-1][2:]
        elif linkindex[i][4][2]==1:
            unitigtoattach=reversecomplement(unitigs[linkindex[i][4][0]-1][2:])
        addcontig(arr1,linkindex[i], unitigtoattach, arr4)
        if contigindex==1:
            break                            
                            
  
def addcontig(arr1, arr2, arr3, arr4):
    temp=[]
    arrs=copy.deepcopy(arr1)
    arrs2=copy.deepcopy(arr4)
    arrs.append([arr2[1], arr1[-1][1]+arr2[2],0,arr2[4][:2]])
    for i in range(0,len(arr3)):
        arrs.append([arr3[i][1], arrs[-1][1]+arr3[i][2],0,arr2[4][:2]])
    for i in range(0,len(arr4)):
        for j in range(len(arr1), len(arrs)):
            if (arrs[arr4[i][0]][0][0]==arrs[j][0][0]+1 and arrs[arr4[i][0]][0][0]%2==0) or (arrs[arr4[i][0]][0][0]==arrs[j][0][0]-1 and arrs[j][0][0]%2==0) :
                if arrs[arr4[i][0]][0][1]==0 and arrs[j][0][1]==1:#should start from 5 end to 3 end               
                    if 3100>=arrs[j][1]-arrs[arr4[i][0]][1]>=1900:
                        if arrs[arr4[i][0]][2]==0 and arrs[j][2]==0:
                            arrs[arr4[i][0]][2]=1
                            arrs[j][2]=1
                            temp.append(i)
                            break
                        elif arrs[arr4[i][0]][0][2]!=arrs[j][2]:
                            print "Something needs to be checked."                            
                    
    for i in range(0, len(temp)):
        arrs2.pop(temp[len(temp)-1-i])
    
    for i in range(len(arr1), len(arrs)):
        for j in range(i+1, len(arrs)):
            if (arrs[i][0][0]==arrs[j][0][0]+1 and arrs[i][0][0]%2==0) or (arrs[i][0][0]==arrs[j][0][0]-1 and arrs[j][0][0]%2==0) :
                if arrs[i][0][1]==0 and arrs[j][0][1]==1:#should start from 5 end to 3 end
                    if 3100>=arrs[j][1]-arrs[i][1]>=1900:
                        if arrs[i][2]==0 and arrs[j][2]==0:
                            arrs[i][2]=1
                            arrs[j][2]=1
                            break
                        elif arrs[i][2]!=arrs[j][2]:
                            print "Something needs to be checked."                    
           
    for i in range(len(arr1), len(arrs)):
        if arrs[i][2]==0:
            arrs2.append([i,arrs[i][1]])     
            
    for i in range(0,len(arrs2)):
        if arrs[len(arrs)-1][1]-arrs2[i][1]>=3100:#make-pair check
            if arrs[arrs2[i][0]-1][2]==1:
                for j in range(arrs2[i][0],len(arrs)):
                    if arrs[j][2]==1:
                        break
                if arrs[j][1]-arrs[arrs2[i][0]-1][1]>500:
                    return                
            if arrs[arrs2[i][0]][3][0]!=3:# Assuming 3 is the only repeat unitig
                return 
        else:                
            break
    contigiteration(arrs,arrs2)



#contig=[a,b,c]        
#a is a read
#b is the read position
#c=1 is if the read satisfies the condition for contig: pairs are within 2400 to 3600
       
checklist_temp=[]

def contigcheckfromzero(arr):
    for i in range(0, len(arr)):
        for j in range(i+1, len(arr)):
            if (arr[i][0][0]==arr[j][0][0]+1 and arr[i][0][0]%2==0) or (arr[i][0][0]==arr[j][0][0]-1 and arr[j][0][0]%2==0):
                if arr[i][0][1]==0 and arr[j][0][1]==1:#should start from 5 end to 3 end
                    if 3100>=arr[j][1]-arr[i][1]>=1900: 
                        if arr[i][2]==0 and arr[j][2]==0:
                            arr[i][2]=1
                            arr[j][2]=1
                            break
                        elif arr[i][2]!=arr[j][2]:
                            print "Something needs to be checked."
    for i in range(0, len(arr)):
        if contig_temp[i][2]==0:
            checklist_temp.append([i,arr[i][1]])
                
                
if len(possiblestartingpoints)>0:
    startingpoint=possiblestartingpoints[0]
    contig_temp=[]
    if startingpoint[1]==0:
        contig_temp.append([unitigs[startingpoint[0]-1][2][0],0,0,startingpoint])
        for i in range(1, len(unitigs[startingpoint[0]-1][2:])):
            contig_temp.append([unitigs[startingpoint[0]-1][2+i][0], contig_temp[i-1][1]+unitigs[startingpoint[0]-1][1+i][2],0,startingpoint])
        contig_temp.append([unitigs[startingpoint[0]-1][2+i][1], contig_temp[i][1]+unitigs[startingpoint[0]-1][2+i][2],0,startingpoint])
        contigcheckfromzero(contig_temp)
        
        contigiteration(contig_temp,checklist_temp)  
    else:
        print "Not implemented yet."
else:
    print "No possible starting points. Think of a middle point."
 
overlap_file = open("lab01.contiginfo", "w")
for i in range(0,len(contigprimitive)):
        print >> overlap_file, contigprimitive[i]
for i in range(0,len(checklist)):
    print >> overlap_file, checklist[i]
overlap_file.close()  

#find mate-pair
contig_pre=[0] * 300
for i in range(0, len(contigprimitive)):
    for j in range(i+1, len(contigprimitive)):
        if (contigprimitive[i][0][0]==contigprimitive[j][0][0]+1 and contigprimitive[i][0][0]%2==0) or (contigprimitive[i][0][0]==contigprimitive[j][0][0]-1 and contigprimitive[j][0][0]%2==0):
            if contigprimitive[i][0][1]==0 and contigprimitive[j][0][1]==1:#should start from 5 end to 3 end
                if 3100>=contigprimitive[j][1]-contigprimitive[i][1]>=1900: 
                    if contig_pre[contigprimitive[i][0][0]-1]==0 and contig_pre[contigprimitive[j][0][0]-1]==0:
                        contig_pre[contigprimitive[i][0][0]-1]=[contigprimitive[i][0], contigprimitive[i][1],contigprimitive[i][1]+499, contigprimitive[i][3]]
                        contig_pre[contigprimitive[j][0][0]-1]=[contigprimitive[j][0], contigprimitive[j][1],contigprimitive[j][1]+499,  contigprimitive[j][3]]
                        break

for i in range(0, len(contig_pre)):
    if contig_pre[i]==0:
        print "error"

overlap_file = open("lab01.extra", "w")
for i in range(0,len(contig_pre)):
        print >> overlap_file, format(contig_pre[i][0][0],'03d'), " ",contig_pre[i][1+contig_pre[i][0][1]]," ", contig_pre[i][2-contig_pre[i][0][1]]
overlap_file.close()

print >> sys.stderr, "Ends at", str(datetime.now())

                    
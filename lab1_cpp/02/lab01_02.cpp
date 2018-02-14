//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

 

#if SAMPLE
const int num_of_reads = 10;
const int read_len = 250;
#else
const int num_of_reads = 300;
const int read_len = 500;
#endif

int read_from_olaps(int list_of_olaps[][5], int location[num_of_reads], vector<vector<int> >  &list_of_exact_olaps){

  char str0[5],str1[5],str2[5],str3[5];
   int num_of_olaps=0;
   int temp0,temp,FB;
   vector<int> temp1(3);
   location[0]=0;

#if SAMPLE
  ifstream infile("sample.olaps");
#else
  ifstream infile("lab01.olaps");
#endif

  if (infile.is_open()) {
    while (infile >> str0>>str1>>str2>>str3){
      temp0=atoi(str3);
      if(temp0==0){
	if (str2[0]=='F')
	  FB=1;
	else
	  FB=0;	
	temp1[0]=atoi(str0)-1;
	temp1[1]=atoi(str1)-1;
	temp1[2]=FB;
	list_of_exact_olaps.push_back(temp1);
      }

      else{
	temp=atoi(str0);
	list_of_olaps[num_of_olaps][0]=temp-1;
	location[temp]=num_of_olaps;

	list_of_olaps[num_of_olaps][1]=atoi(str1)-1;
      
	if (str2[0]=='F')
	  list_of_olaps[num_of_olaps][2]=1;
	else
	  list_of_olaps[num_of_olaps][2]=0;

	list_of_olaps[num_of_olaps][3]=atoi(str3);
	list_of_olaps[num_of_olaps][4]=0;

	num_of_olaps++;
      }
      
    }
  }
  infile.close();
 
  for (int i=1;i<num_of_reads;i++){
    if (location[i]==0)
      location[i]=location[i-1];
    else location[i]+=1;
  }
  
  return num_of_olaps;
}

int to_delete(int i, int j, int k, int l, int list_of_olaps[][5]){

  
  return -1;
}

void delete_edges(int list_of_olaps[][5],int num_of_olaps, int location[num_of_reads]){

  /*i the number of read1, j the location of read1, read2 pair, k the location of read1, read2 pair, l, the number of
		read2, read2 pair*/

  int i,j,k,l;
  int deleted=1;
  int todelete;
  
  for (i=0;i<num_of_reads-2;i++){
    for (j=location[i];j<location[i+1]-1;j++){
      for (k=j+1;k<location[i+1];k++){
	for (l=location[list_of_olaps[j][1]];l<location[list_of_olaps[j][1]+1];l++){
	  if(list_of_olaps[l][1]==list_of_olaps[k][1]){
	    todelete=to_delete(i,j,k,l,list_of_olaps);
	    list_of_olaps[todelete][4]=1;
	  }
	  else if (list_of_olaps[l][1]>list_of_olaps[k][1])
	    break;
	  i, list_of_olaps[j][1] list_of_olaps[k][1] list_of_olaps[l][0] list_of_olaps[1]
	}
      }
    }
  }
}

int main(){

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][5];//read1,read2,F/R,length,deleted

  vector<vector<int> > list_of_exact_olaps(0, vector<int>(3));//save exact overlaps
  
  int location[num_of_reads];//starting point of the read 1 in the list_of_olaps

  int num_of_olaps=read_from_olaps(list_of_olaps,location,list_of_exact_olaps);

  delete_edges(list_of_olaps,num_of_olaps,location);
  
  return 0;
}

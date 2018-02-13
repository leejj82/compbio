//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;
typedef vector<int> int_vec;
 

#if SAMPLE
const int num_of_reads = 10;
const int read_len = 250;
#else
const int num_of_reads = 300;
const int read_len = 500;
#endif

int read_from_olaps(int list_of_olaps[][5], int location[num_of_reads]){

  char str0[5],str1[5],str2[5],str3[5];
   int num_of_olaps=0;
   int temp;
   location[0]=0;

#if SAMPLE
  ifstream infile("sample.olaps");
#else
  ifstream infile("lab01.olaps");
#endif

  if (infile.is_open()) {
    while (infile >> str0>>str1>>str2>>str3){
      
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
  infile.close();
 
  for (int i=1;i<num_of_reads;i++){
    if (location[i]==0)
      location[i]=location[i-1];
    else location[i]+=1;
  }
  
  return num_of_olaps;
}

void delete_edges(int list_of_olaps[][5],int num_of_olaps, int location[num_of_reads]){
  int i,j,k;
  for (i=0;i<num_of_olaps-2;i++){
    for (j=i+1;j<i+num_of_reads;j++){
    }
  }
}

int main(){

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][5];//read1,read2,F/R,length,deleted
  int location[num_of_reads];//starting point of the read 1 in the list_of_olaps

  int num_of_olaps= num_of_olaps=read_from_olaps(list_of_olaps,location);

  delete_edges(list_of_olaps,num_of_olaps,location);

  int_vec v;
  cout<<v.size()<<v.capacity();
  
  
  return 0;
}

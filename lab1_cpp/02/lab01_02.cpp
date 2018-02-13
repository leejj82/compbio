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

int read_from_olaps(int list_of_olaps[][5]){

  char str0[5],str1[5],str2[5],str3[5];
   int num_of_olaps=0;

#if SAMPLE
  ifstream infile("sample.olaps");
#else
  ifstream infile("lab01.olaps");
#endif

  if (infile.is_open()) {
    while (infile >> str0>>str1>>str2>>str3){
	list_of_olaps[num_of_olaps][0]=atoi(str0)-1;
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
  return num_of_olaps;
}

int main(){

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][5];
  int num_of_olaps;
  
  num_of_olaps=read_from_olaps(list_of_olaps);

  cout<<num_of_olaps<<"\n";
  for (int j=0;j<100;j++){
  for (int i=0;i<5;i++)
    cout<<list_of_olaps[j][i]<<" ";
  cout<<"\n";}

  
  int_vec v;
  cout<<v.size()<<v.capacity();
  
  
  return 0;
}

//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <iomanip>
#include <assert.h>

using namespace std;

 

#if SAMPLE
const int num_of_reads = 10;
const int read_len = 250;
#else
const int num_of_reads = 300;
const int read_len = 500;
#endif

int read_from_olaps(int list_of_olaps[][6], int location[num_of_reads], vector<vector<int> >  &list_of_exact_olaps, int &num_of_exact_olaps){

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
	num_of_exact_olaps+=1;
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

	list_of_olaps[num_of_olaps][3]=temp0;
	if (temp0 >0)
	  list_of_olaps[num_of_olaps][4]=1;
	else
	  list_of_olaps[num_of_olaps][4]=0;
	list_of_olaps[num_of_olaps][5]=0;
	
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

int record_edge_to_delete(int list_of_olaps[][6],int num_of_olaps, int location[num_of_reads]){

  /*i the number of read1, j the location of read1, read2 pair, k the location of read1, read2 pair, l, the number of
		read2, read2 pair*/

  int i,j,k,l;
  int deleted=1;
  int arrow0,arrow1,arrow2;
  int count=0;
  
  for (i=0;i<num_of_reads-2;i++){
    for (j=location[i];j<location[i+1]-1;j++){
      for (k=j+1;k<location[i+1];k++){
	for (l=location[list_of_olaps[j][1]];l<location[list_of_olaps[j][1]+1];l++){
      	  if(list_of_olaps[k][1]==list_of_olaps[l][1]){
    	    if((list_of_olaps[j][2]+list_of_olaps[k][2]+list_of_olaps[l][2])%2==1){
	      arrow0=list_of_olaps[j][4];
	      arrow1=list_of_olaps[k][4];
	      arrow2=list_of_olaps[l][4];
	      if (list_of_olaps[j][2]==0)
		arrow2=1-arrow2;
	      if (arrow0==arrow1 && arrow1==arrow2 && list_of_olaps[k][5]==0){
		list_of_olaps[k][5]=1;
		count++;
	      }
	      else if (arrow0==arrow1 && arrow0!=arrow2 && list_of_olaps[j][5]==0){
		list_of_olaps[j][5]=1;
		count++;
	      }
	      else if (arrow0!=arrow1 && arrow1==arrow2 && list_of_olaps[l][5]==0){
		list_of_olaps[l][5]=1;
		count++;
	      }
	    }
	  }
	  else if (list_of_olaps[l][1]>list_of_olaps[k][1])
	    break;
	}
      }
    }
  }
  return count;
}

void find_unitigs(int location[num_of_reads],int num_of_olaps,int list_of_olaps[][6],vector<vector<vector<int> > > &edges_for_nodes){

  int i,j=0,k;
  int count[2];
  vector<int> temp(2);
  vector<vector<int> > temp_node(2,vector<int>(0));
  
  for (i=0;i<num_of_reads;i++){
    count[0]=0;
    count[1]=0;
    for (j=location[i];j<location[i+1];j++){
      if (list_of_olaps[j][5]==0){
	if (list_of_olaps[j][4]==1){
     	  temp_node[0].push_back(list_of_olaps[j][0]);
	  temp_node[0].push_back(list_of_olaps[j][1]);	  
	  count[0]++;
	}
	else {
     	  temp_node[1].push_back(list_of_olaps[j][0]);
	  temp_node[1].push_back(list_of_olaps[j][1]);	  
	  count[1]++;
	}
      }
    }
    edges_for_nodes.push_back(temp_node);
    temp_node[0].resize(0);
    temp_node[1].resize(0);
  }
}   




int main(){

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][6];//read1,read2,F/R(1/0),olap_length,deleted,direction of arrow

  vector<vector<int> > list_of_exact_olaps;//save exact overlaps
  int num_of_exact_olaps=0;
  
  int location[num_of_reads];//starting point of the read 1 in the list_of_olaps

  int num_of_olaps=read_from_olaps(list_of_olaps,location,list_of_exact_olaps,num_of_exact_olaps);

  int num_of_edges_to_delete= record_edge_to_delete(list_of_olaps,num_of_olaps,location);

  int num_of_edges_for_unitigs=num_of_olaps-num_of_edges_to_delete;

  vector<vector<vector<int> > > edges_for_nodes(num_of_reads);
  edges_for_nodes.clear();
  
  find_unitigs(location, num_of_olaps,list_of_olaps,edges_for_nodes);
  cout<<num_of_edges_to_delete<<" "<<num_of_edges_for_unitigs<<" ";

  cout<<list_of_exact_olaps.size()<<" ";
  for (int i=0;i<list_of_exact_olaps.size();i++){
    cout<<list_of_exact_olaps[i].size()<<" ";
    for (int j=0;j<list_of_exact_olaps[i].size();j++){
      cout<<list_of_exact_olaps[i][j]<<" ";
    }
  }

  // cout << setw(10) << 77 << endl;

  int i,j,k;

  FILE * pFile;
  pFile = fopen("lab01.list_edges", "w");
  assert(edges_for_nodes.size() == num_of_reads);  
  for(i=0; i < edges_for_nodes.size(); i++) {
    assert(edges_for_nodes[i].size() == 2);
    for(j=0; j < edges_for_nodes[i].size(); j++) {
      // cout<< edges_for_nodes[i][j].size() <<"\n";
      /*      for (k=0;k<edges_for_nodes[i][j].size();k++){
	      fprintf (pFile, "%d\n",edges_for_nodes[i][j][k]);
      }*/
    }
  }
 

  fclose (pFile);
  
  return 0;
}

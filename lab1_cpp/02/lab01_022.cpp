//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
<<<<<<< HEAD
#include<list>
=======
#include <iomanip>
#include <assert.h>
>>>>>>> 6a9ac841bfa03422c51d928f501cc7e70a5c1dd2

using namespace std;

 

#if SAMPLE
const int num_of_reads = 10;
const int read_len = 250;
#else
const int num_of_reads = 300;
const int read_len = 500;
#endif

int read_from_olaps(int list_of_olaps[][6], int location[num_of_reads], vector<vector<int> >  &list_of_exact_olaps, int &num_of_exact_olaps,int &num_of_edges_to_delete){

  char str0[5],str1[5],str2[5],str3[5];
  int i,j,temp0,temp,FB,num_of_olaps=0;
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
 
  for (i=1;i<num_of_reads;i++){
    if (location[i]==0)
      location[i]=location[i-1];
    else location[i]+=1;
  }

  for (i=0;i<num_of_exact_olaps;i++){
    for (j=0;j<num_of_olaps;j++){
      if(list_of_exact_olaps[i][1]==list_of_olaps[j][0] || list_of_exact_olaps[i][1]==list_of_olaps[j][1]){
	list_of_olaps[j][5]=1;
	num_of_edges_to_delete++;
      }
    }
  }

  
  return num_of_olaps;
}

void record_edge_to_delete(int list_of_olaps[][6],int num_of_olaps, int location[num_of_reads],int &num_of_edges_to_delete){

  /*i the number of read1, j the location of read1, read2 pair, k the location of read1, read2 pair, l, the number of
		read2, read2 pair*/

  int i,j,k,l;
  int deleted=1;
  int arrow0,arrow1,arrow2;
   
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
		list_of_olaps[k][5]=deleted;
		num_of_edges_to_delete++;
	      }
	      else if (arrow0==arrow1 && arrow0!=arrow2 && list_of_olaps[j][5]==0){
		list_of_olaps[j][5]=deleted;
		num_of_edges_to_delete++;
	      }
	      else if (arrow0!=arrow1 && arrow1==arrow2 && list_of_olaps[l][5]==0){
		list_of_olaps[l][5]=deleted;
		num_of_edges_to_delete++;
      	      }
	    }
	  }
	  else if (list_of_olaps[l][1]>list_of_olaps[k][1])
	    break;
	}
      }
    }
  }
}

void set_up_viable_edges(int location[num_of_reads], int list_of_olaps[][6], vector<vector<vector<int> > > &edges_for_nodes,int edges_for_nodes_index[][4],  vector<vector<int> > &list_of_exact_olaps, int num_of_exact_olaps){

  int i,j;
  
  for (i=0;i<num_of_reads;i++){
    for (j=location[i];j<location[i+1];j++){
      if (list_of_olaps[j][5]==0){
	if (list_of_olaps[j][4]==1){
	  edges_for_nodes[i][0].push_back(list_of_olaps[j][1]);
	  edges_for_nodes[i][0].push_back(list_of_olaps[j][2]);
	  edges_for_nodes[i][0].push_back(list_of_olaps[j][3]);
	  edges_for_nodes[i][0].push_back(0);
	  edges_for_nodes_index[i][0]+=1;
	  edges_for_nodes_index[i][1]+=1;

	  if (list_of_olaps[j][2]==1){
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][2]);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(-list_of_olaps[j][3]);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(0);
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][2]+=1;
	  }
	  else{
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][2]);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][3]);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(0);
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][1]+=1;
	    
	  }
	}
	else {
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][1]);
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][2]);
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][3]);
	  edges_for_nodes[i][1].push_back(0);
	  edges_for_nodes_index[i][0]+=1;
	  edges_for_nodes_index[i][2]+=1;
	  

	  if (list_of_olaps[j][2]==1){
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][2]);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(-list_of_olaps[j][3]);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(0);
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][1]+=1;
	  }
	  else{
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][2]);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][3]);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(0);
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][2]+=1;
	  }
	}
      }
    }
  }

  for (i=0;i<num_of_exact_olaps;i++){
    edges_for_nodes_index[list_of_exact_olaps[i][1]][3]+=1;
    cout<<list_of_exact_olaps[i][1];
  }
}

/*
void find_a_unitig(vector<list<vector<int> > > &unitigs, int u_count, vector<vector<vector<int> > > &edges_for_nodes, int &index,int edges_for_nodes_index[][4]){

  int i;
  for (i=index;i<num_of_reads;i++){
    if (edges_for_nodes_index[i][3]==0){//not used
      if (edges_for_nodes_index[i][1]==1){
	unitigs[u_count].push_back(edges_for_nodes[i][0]);
      }
      if (edges_for_nodes_index[i][2]==1 {
	
	}

	

	}
    }
  } 
}
*/


void find_unitigs(vector<list<vector<int> > > &unitigs, vector<vector<vector<int> > > &edges_for_nodes, int edges_for_nodes_index[][4]){
  /*
  int i,u_count=0, index=0;
  int front=0, back=1
  while(index<num_of_reads){
    for (i=index;i<num_of_reads;i++){
      if (edges_for_nodes_index[i][3]==0){//not used
	
	list<vector<int> > temp;
	unitigs.push_back(temp);

	if (edges_for_nodes_index[i][1]==1 && edges_for_nodes_index[i][2]==1){
	  insert_unitig(front, unitigs, u_count, edges_for_nodes[i][1]);
	  insert_unitig(back, unitigs, u_count, edges_for_nodes[i][0]);
	  
	  unitigs[u_count].push_front(edges_for_nodes[i][1]);
	  unitigs[u_count].push_back(edges_for_nodes[i][0]);
	 
	  //	  find_a_unitig(unitigs,u_count,edges_for_nodes,index, edges_for_nodes_index);
	}
	else{
	  if (edges_for_nodes_index[i][1]==1){

	  }
	  else if (edges_for_nodes_index[i][2]==1){

	  }
	  else {

	  }
	}
      }
    }
  }*/
}   




int main(){

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][6];//read1,read2,F/R(1/0),olap_length,deleted,direction of arrow

  vector<vector<int> > list_of_exact_olaps;//save exact overlaps
  int num_of_exact_olaps=0;
  int num_of_edges_to_delete=0;
  int location[num_of_reads];//starting point of the read 1 in the list_of_olaps

  int num_of_olaps=read_from_olaps(list_of_olaps,location,list_of_exact_olaps,num_of_exact_olaps,num_of_edges_to_delete);

  record_edge_to_delete(list_of_olaps,num_of_olaps,location, num_of_edges_to_delete);

  int num_of_edges_for_unitigs=num_of_olaps-num_of_edges_to_delete;

<<<<<<< HEAD
  vector<vector<vector<int> > > edges_for_nodes(num_of_reads, vector<vector<int> >(2));//for each node, save outgoing edges, incoming edges

  int edges_for_nodes_index[num_of_reads][4]={{0}}; //save #of total edges, # of outgoing edges,  # of incoming edges, #used or not
  
  set_up_viable_edges(location, list_of_olaps, edges_for_nodes,edges_for_nodes_index, list_of_exact_olaps, num_of_exact_olaps);

  vector<list<vector<int> > > unitigs;

  find_unitigs(unitigs,edges_for_nodes,edges_for_nodes_index);
=======
  vector<vector<vector<int> > > edges_for_nodes(num_of_reads);
  edges_for_nodes.clear();
  
  find_unitigs(location, num_of_olaps,list_of_olaps,edges_for_nodes);
  cout<<num_of_edges_to_delete<<" "<<num_of_edges_for_unitigs<<" ";
>>>>>>> 6a9ac841bfa03422c51d928f501cc7e70a5c1dd2







  

  cout<<num_of_edges_to_delete<<" "<<num_of_edges_for_unitigs<<" ";
  cout<<list_of_exact_olaps.size()<<" ";
  for (int i=0;i<list_of_exact_olaps.size();i++){
    cout<<list_of_exact_olaps[i].size()<<" ";
    for (int j=0;j<list_of_exact_olaps[i].size();j++){
      cout<<list_of_exact_olaps[i][j]<<" ";
    }
  }
<<<<<<< HEAD
  cout<<endl;
=======

  // cout << setw(10) << 77 << endl;
>>>>>>> 6a9ac841bfa03422c51d928f501cc7e70a5c1dd2

  int i,j,k;

  FILE * pFile;
<<<<<<< HEAD
  pFile = fopen ("lab01.list_edges","w");

  for (i=0;i<num_of_reads;i++){
    for (j=0;j<2;j++){
      fprintf (pFile, "%d  ",i);
	  
      fprintf (pFile, "%d  ",edges_for_nodes_index[i][0]);
      fprintf (pFile, "%d  ",edges_for_nodes_index[i][j+1]);
      fprintf (pFile, "%d  ",edges_for_nodes_index[i][3]);
      
      
      for (k=0;k<edges_for_nodes[i][j].size();k++){
	fprintf (pFile, "%d  ",edges_for_nodes[i][j][k]);
      }

      fprintf (pFile, "\n");
=======
  pFile = fopen("lab01.list_edges", "w");
  assert(edges_for_nodes.size() == num_of_reads);  
  for(i=0; i < edges_for_nodes.size(); i++) {
    assert(edges_for_nodes[i].size() == 2);
    for(j=0; j < edges_for_nodes[i].size(); j++) {
      // cout<< edges_for_nodes[i][j].size() <<"\n";
      /*      for (k=0;k<edges_for_nodes[i][j].size();k++){
	      fprintf (pFile, "%d\n",edges_for_nodes[i][j][k]);
      }*/
>>>>>>> 6a9ac841bfa03422c51d928f501cc7e70a5c1dd2
    }
  }
 

  fclose (pFile);
  
  return 0;
}

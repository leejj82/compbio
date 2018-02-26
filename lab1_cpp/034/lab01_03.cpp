//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<list>

using namespace std;

 
#if SAMPLE
const int num_of_reads = 10;
const int read_len = 250;
const int l_bd_mp= 950;
const int u_bd_mp= 1000;
#else
const int num_of_reads = 300;
const int read_len = 500;
const int l_bd_mp= 1900;
const int u_bd_mp= 3100;
#endif

void read_from_fasta(char list_of_reads[][read_len + 1]){

  char str[500];
   int i=0;

#if SAMPLE
  ifstream infile("sample.fasta");
#else
  ifstream infile("lab01.fasta");
#endif

  if (infile.is_open()) {
    while (infile >> str)
      {
	if (str[0]=='>')
	 {
	   i=atoi(str+1)-1;
	 }
	else
	  {
	    strcat(list_of_reads[i],str);
	    }
      }
  }
  infile.close();
}

void reversecomplement(char list[][read_len + 1], char list_rc[][read_len + 1]){
  for (int j=0;j<num_of_reads;j++)
    {
      for (int i=0;i<read_len;i++)
	{
	  if (list[j][read_len-i-1]=='a')
	    {
	      list_rc[j][i]='t';
	    }
	  else if (list[j][read_len-i-1]=='t')
	    {
	      list_rc[j][i]='a'; 
	    }
	  else if (list[j][read_len-i-1]=='c')
	    {
	      list_rc[j][i]='g';
	    }
	  else if (list[j][read_len-i-1]=='g')
	    {
	      list_rc[j][i]='c';
	    }
	}
      list_rc[j][read_len]='\0';
    }
 }

int read_from_olaps(int list_of_olaps[][6], int location[num_of_reads+1], vector<vector<int> >  &list_of_exact_olaps, int &num_of_exact_olaps,int &num_of_edges_to_delete){

  char str0[5],str1[5],str2[5],str3[5];
  int i,j,temp0,temp,FB,num_of_olaps=0;
  vector<int> temp1(3);


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
  location[num_of_reads]=num_of_olaps;
  
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

void record_edge_to_delete(int list_of_olaps[][6],int num_of_olaps, int location[num_of_reads+1],int &num_of_edges_to_delete){

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
  
  for (i=0;i<num_of_reads;i++){//read1,read2,F/R(1/0),olap_length,location,deleted
    for (j=location[i];j<location[i+1];j++){
      if (list_of_olaps[j][5]==0){ //not deleted
	if (list_of_olaps[j][4]==1){ //location of first read front
	  edges_for_nodes[i][1].push_back(i);
	  edges_for_nodes[i][1].push_back(1);
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][1]);
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][2]);
	  edges_for_nodes[i][1].push_back(list_of_olaps[j][3]);
	
	  edges_for_nodes_index[i][0]+=1;
	  edges_for_nodes_index[i][2]+=1;

	  if (list_of_olaps[j][2]==1){//forward
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(1);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][1]);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(1);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][3]);
	    
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][1]+=1;
	  }
	  else{//reverse complement
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][1]);
     	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(1);
     	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(0);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][3]);
	    
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][2]+=1;
	  }
	}
	else { //location of first read back
	  edges_for_nodes[i][0].push_back(list_of_olaps[j][1]);
	  edges_for_nodes[i][0].push_back(list_of_olaps[j][2]);
	  edges_for_nodes[i][0].push_back(i);
	  edges_for_nodes[i][0].push_back(1);
	  edges_for_nodes[i][0].push_back(-list_of_olaps[j][3]);
	
	  edges_for_nodes_index[i][0]+=1;
	  edges_for_nodes_index[i][1]+=1;

	  if (list_of_olaps[j][2]==1){//forward
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(list_of_olaps[j][1]);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(1);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(1);
	    edges_for_nodes[list_of_olaps[j][1]][1].push_back(-list_of_olaps[j][3]);
	    
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][2]+=1;
	  }
	  else{//reverse complement
  	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(i);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(0);
	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(list_of_olaps[j][1]);
     	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(1);
   	    edges_for_nodes[list_of_olaps[j][1]][0].push_back(-list_of_olaps[j][3]);
	    
	    edges_for_nodes_index[list_of_olaps[j][1]][0]+=1;
	    edges_for_nodes_index[list_of_olaps[j][1]][1]+=1;
	  }
	}
      }
    }
  }

  for (i=0;i<num_of_exact_olaps;i++){
    edges_for_nodes_index[list_of_exact_olaps[i][1]][3]+=1;
  }
}
 
void set_up_edges_RC(vector<vector<vector<int> > > &edges_for_nodes, vector<vector<vector<int> > > &edges_for_nodes_RC){

  int temp;
  for (int i=0;i<edges_for_nodes.size();i++){
    for (int j=0;j<2;j++){
      if (edges_for_nodes[i][j].size()>0){
	edges_for_nodes_RC[i][1-j]=edges_for_nodes[i][j];
	for (int k=0;k<edges_for_nodes[i][j].size()/5;k++){
	  temp=edges_for_nodes_RC[i][1-j][5*k];
	  edges_for_nodes_RC[i][1-j][5*k]=edges_for_nodes_RC[i][1-j][5*k+2];
	  edges_for_nodes_RC[i][1-j][5*k+2]=temp;

	  temp=edges_for_nodes_RC[i][1-j][5*k+1];    
	  edges_for_nodes_RC[i][1-j][5*k+1]=1-edges_for_nodes_RC[i][1-j][5*k+3];
	  edges_for_nodes_RC[i][1-j][5*k+3]=1-temp;
	}
      }
    }
  }
}


void previous_read(int starting_point, int FB,  vector<vector<vector<int> > > &edges_for_nodes,  vector<vector<vector<int> > > &edges_for_nodes_RC, int edges_for_nodes_index[][4], vector<vector<int> > &unitig_front){

  int i,used=1;
  int Forward=1, RC=0;
  int previous_node;
  int array[5]={-1,-1,-1,-1,-1};
  vector<int> null_vector(array,array+5);
  
  if (FB==Forward){ //current node is in the forward order
    if (edges_for_nodes_index[starting_point][1]>=1){  

      unitig_front.push_back(edges_for_nodes[starting_point][0]);

      if (edges_for_nodes_index[starting_point][1]==1){//there is exactly one edge connecting with previous node
	if (edges_for_nodes[starting_point][0][1]==Forward){//the previous node is in the forward order

	  previous_node=edges_for_nodes[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][2]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }      
	}
	else{//the previous node is in the reverse complement order
      
	  previous_node=edges_for_nodes[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][1]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }
	}
      }
    }
    else { //there is no edge to the previous node=>attach null vector
      unitig_front.push_back(null_vector);
    }
  }
  else{ //current node is in the reverse complement order
    
    if (edges_for_nodes_index[starting_point][2]>=1){  

      unitig_front.push_back(edges_for_nodes_RC[starting_point][0]);

      if (edges_for_nodes_index[starting_point][2]==1){//there is exactly one edge connecting with previous node
	if (edges_for_nodes_RC[starting_point][0][1]==Forward){//the previous node is in the forward order

	  previous_node=edges_for_nodes_RC[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][2]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }      
	}
	else{//the previous node is in the reverse complement order
      
	  previous_node=edges_for_nodes_RC[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][1]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }
	}
      }
    }
    else { //there is no edge to the next node=>attach null vector
      unitig_front.push_back(null_vector);
    }
  }
}

void next_read(int starting_point, int FB,  vector<vector<vector<int> > > &edges_for_nodes,  vector<vector<vector<int> > > &edges_for_nodes_RC, int edges_for_nodes_index[][4], vector<vector<int> > &unitig_back){
  
  int i,used=1;
  int Forward=1, RC=0;
  int next_node;

  int array[5]={-1,-1,-1,-1,-1};
  vector<int> null_vector(array,array+5);
    
  if (FB==Forward){ //current node is in the forward order
    if (edges_for_nodes_index[starting_point][2]>=1){//next node

      unitig_back.push_back(edges_for_nodes[starting_point][1]);

      if (edges_for_nodes_index[starting_point][2]==1){//there is exactly one edge connecting with next node
	if (edges_for_nodes[starting_point][1][3]==Forward){//the next node is in the forward order

	  next_node=edges_for_nodes[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][1]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
	else{//the next node is in the reverse complement order
      
	  next_node=edges_for_nodes[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][2]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
      }
    }
    else { //there is no edge to the previous node=>attach null vector
      unitig_back.push_back(null_vector);
    }
  }
  else{ //current node is in the reverse complement order
    if (edges_for_nodes_index[starting_point][1]>=1){//next node

      unitig_back.push_back(edges_for_nodes_RC[starting_point][1]);

      if (edges_for_nodes_index[starting_point][1]==1){//there is exactly one edge connecting with next node
	if (edges_for_nodes_RC[starting_point][1][3]==Forward){//the next node is in the forward order

	  next_node=edges_for_nodes_RC[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][1]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
	else{//the next node is in the reverse complement order
      
	  next_node=edges_for_nodes_RC[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][2]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
      }
    }
    else { //there is no edge to the next node=>attach null vector
      unitig_back.push_back(null_vector);
    }
  }
}

void find_a_unitig(int &starting_point, vector<vector<vector<int> > > &edges_for_nodes, vector<vector<vector<int> > > &edges_for_nodes_RC, int edges_for_nodes_index[][4], vector<vector<vector<int> > >  &unitigs,vector<vector<int> > &unitigs_info, vector<vector<int> > &list_of_exact_olaps, int num_of_exact_olaps, int &read_exact_match_count){

  int i,j,used=1;
  vector<vector<int> > unitig_front;
  vector<vector<int> > unitig_back;
  vector<vector<int> > unitig;
  int Forward=1, RC=0;
  int previous_node, next_node;

  int array[5]={-1,-1,-1,-1,-1};
  vector<int> null_vector(array,array+5);
  vector<int> temp(2);

  int sum;
  
  if (edges_for_nodes_index[starting_point][0]>0){ //there exists at least one edge connected to the node
  
    if (edges_for_nodes_index[starting_point][1]>=1){  //previous node

      unitig_front.push_back(edges_for_nodes[starting_point][0]);

      if (edges_for_nodes_index[starting_point][1]==1){//there is exactly one edge connecting with previous node
	if (edges_for_nodes[starting_point][0][1]==Forward){//the previous node is in the forward order

	  previous_node=edges_for_nodes[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][2]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }      
	}
	else{//the previous node is in the reverse complement order
      
	  previous_node=edges_for_nodes[starting_point][0][0];

	  if (edges_for_nodes_index[previous_node][1]==1 ){//the previous node has exactly one outgoing edge
	    edges_for_nodes_index[previous_node][3]=used;
	    previous_read(previous_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_front);
	  }
	}
      }
    }
    else { //there is no edge to the previous node=>attach null vector
      unitig_front.push_back(null_vector);
    }

    if (edges_for_nodes_index[starting_point][2]>=1){//next node

      unitig_back.push_back(edges_for_nodes[starting_point][1]);

      if (edges_for_nodes_index[starting_point][2]==1){//there is exactly one edge connecting with next node
	if (edges_for_nodes[starting_point][1][3]==Forward){//the next node is in the forward order

	  next_node=edges_for_nodes[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][1]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, Forward, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
	else{//the next node is in the reverse complement order
      
	  next_node=edges_for_nodes[starting_point][1][2];

	  if (edges_for_nodes_index[next_node][2]==1 ){//the next node has exactly one incoming edge
	    edges_for_nodes_index[next_node][3]=used;
	    next_read(next_node, RC, edges_for_nodes,edges_for_nodes_RC,edges_for_nodes_index,unitig_back);
	  }
	}
      }
    }
    else { //there is no edge to the next node=>attach null vector
      unitig_back.push_back(null_vector);
    }
    
    reverse(unitig_front.begin(),unitig_front.end());//reorganize unitig
    unitig_front.insert(unitig_front.end(),unitig_back.begin(),unitig_back.end());

    if (read_exact_match_count<num_of_exact_olaps){
      for (i=1;i<unitig_front.size();i++){//insert exactly identical reads
	for (j=0;j<num_of_exact_olaps;j++){
	  if(unitig_front[i][0]==list_of_exact_olaps[j][0]){

	    read_exact_match_count++;
	    if (unitig_front[i][0]!=-1){//not the inserted null vector
	      if(unitig_front[i][1]==1){//front order
	   
		vector<int> vec;
		vec.push_back(unitig_front[i][0]);
		vec.push_back(1);
		vec.push_back(list_of_exact_olaps[j][1]);
		vec.push_back(list_of_exact_olaps[j][2]);
		vec.push_back(0);

		unitig_front[i][0]=list_of_exact_olaps[j][1];
		unitig_front[i][1]=list_of_exact_olaps[j][2];
		unitig_front.insert(unitig_front.begin()+i,vec);
	      }
	      else {//reverse complement order
	   
		vector<int> vec;
		vec.push_back(unitig_front[i][0]);
		vec.push_back(0);
		vec.push_back(list_of_exact_olaps[j][1]);
		vec.push_back(1-list_of_exact_olaps[j][2]);
		vec.push_back(0);

		unitig_front[i][0]=list_of_exact_olaps[j][1];
		unitig_front[i][1]=1-list_of_exact_olaps[j][2];
	    
		unitig_front.insert(unitig_front.begin()+i,vec);
	      }
	    }
	  }
	}
      }
    }

    temp[0]=(int)(unitig_front.size())-1; //# of reads in a unitig
    sum=read_len;//sum is the length of each unitig
    for (i=1;i<temp[0];i++){
      sum+=unitig_front[i][4];
    }
    temp[1]=sum;

    unitigs_info.push_back(temp);//record # of reads and length for each unitig
    unitigs.push_back(unitig_front);//input in the set of unitigs


   
  }

  else {//single node case
    cout<<"There exists a single node. Not implemented yet.";
  }  
}
  

int find_unitigs(vector<vector<vector<int> > >  &unitigs, vector<vector<int> > &unitigs_info,  vector<vector<vector<int> > > &edges_for_nodes, vector<vector<vector<int> > > &edges_for_nodes_RC, int edges_for_nodes_index[][4], vector<vector<int> > &list_of_exact_olaps, int num_of_exact_olaps){
  
  int not_used=0,used=1;
  int front=1,back=0;
  int read_exact_match_count=0;
  
  for(int starting_point=0;starting_point<num_of_reads;starting_point++){
    if (edges_for_nodes_index[starting_point][3]==not_used){//node not used
      edges_for_nodes_index[starting_point][3]=used;
      find_a_unitig(starting_point,edges_for_nodes, edges_for_nodes_RC,edges_for_nodes_index, unitigs, unitigs_info, list_of_exact_olaps, num_of_exact_olaps, read_exact_match_count);
    }
  }
}

int count_the_num_of_connections(vector<vector<vector<int> > > &unitigs,int &num_of_unitigs,vector<vector<int> > &unitigs_info, int unitigs_con_count[][3]){

  int i;
  int total_count=0;
  
  for (i=0;i<num_of_unitigs;i++){
    
    if (unitigs[i][0][0]!=-1){
      unitigs_con_count[i][1]=unitigs[i][0].size()/5;
    } 
    else unitigs_con_count[i][1]=0;
    if(unitigs[i][unitigs_info[i][0]][0]!=-1){
      unitigs_con_count[i][2]=unitigs[i][unitigs_info[i][0]].size()/5;
    }
    else unitigs_con_count[i][2]=0;

    unitigs_con_count[i][0]=unitigs_con_count[i][1]+unitigs_con_count[i][2];
    total_count+=unitigs_con_count[i][0];
  }
  return total_count;
}

void setup_unis(vector<vector<vector<int> > > &unitigs,vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info ){

  int i,j,k;
  
  for (i=0;i<num_of_unitigs;i++){
    unis[i].resize(unitigs_info[i][0], vector<int>(3));  
    unis_RC[i].resize(unitigs_info[i][0], vector<int>(3));

    unis[i][0][0]=unitigs[i][1][0];
    unis[i][0][1]=unitigs[i][1][1];
    unis[i][0][2]=0;

    unis_RC[i][0][0]=unitigs[i][unitigs_info[i][0]-1][2];
    unis_RC[i][0][1]=1-unitigs[i][unitigs_info[i][0]-1][3];
    unis_RC[i][0][2]=0;
  
    for (j=1;j<unitigs_info[i][0];j++){     
      unis[i][j][0]=unitigs[i][j][2];
      unis[i][j][1]=unitigs[i][j][3];
      unis[i][j][2]=unis[i][j-1][2]+unitigs[i][j][4];

      unis_RC[i][j][0]=unitigs[i][unitigs_info[i][0]-1-j][2];
      unis_RC[i][j][1]=1-unitigs[i][unitigs_info[i][0]-1-j][3];
      unis_RC[i][j][2]=unis_RC[i][j-1][2]+unitigs[i][unitigs_info[i][0]-j][4];
    }
  }
}
 
void connected_unitigs(vector<vector<vector<int> > > &unitigs,int &num_of_unitigs, vector<vector<int> > &unitigs_info, int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con){

  int i,j,k,l,m,n; 
  vector<int> temp(5);
  int T;

  for (i=0;i<num_of_unitigs-1;i++){
    for (j=0;j<2;j++){
      for (k=0;k<unitigs_con_count[i][j+1];k++){
	for(l=i+1;l<num_of_unitigs;l++){
	  for (m=0;m<2;m++){
	    for (n=0;n<unitigs_con_count[l][m+1];n++){
	     
	      if(equal(unitigs[i][j*unitigs_info[i][0]].begin()+k*5,unitigs[i][j*unitigs_info[i][0]].begin()+(k+1)*5,unitigs[l][m*unitigs_info[l][0]].begin()+n*5) && m!=j && unitigs[i][j*unitigs_info[i][0]][k*5]!=-1 ){
		temp[0]=(j*i+(1-j)*l);
		temp[1]=1;//Forward=1, Reverse Complement=0
		temp[2]=(j*l+(1-j)*i);
		temp[3]=1;;//Forward=1, Reverse Complement=0
		temp[4]=unitigs[i][j*unitigs_info[i][0]][4+k*5];
						  
		unitigs_con[i][j].push_back(temp);
		unitigs_con[l][m].push_back(temp);
	      }
	      else if (unitigs[i][j*unitigs_info[i][0]][k*5]==unitigs[l][m*unitigs_info[l][0]][2+n*5] && unitigs[i][j*unitigs_info[i][0]][1+k*5]==1-unitigs[l][m*unitigs_info[l][0]][3+n*5] && unitigs[i][j*unitigs_info[i][0]][2+k*5]==unitigs[l][m*unitigs_info[l][0]][0+n*5]  && unitigs[i][j*unitigs_info[i][0]][3+k*5]==1-unitigs[l][m*unitigs_info[l][0]][1+n*5] && unitigs[i][j*unitigs_info[i][0]][4+k*5]==unitigs[l][m*unitigs_info[l][0]][4+n*5] && m==j && unitigs[i][j*unitigs_info[i][0]][k*5]!=-1){

		temp[0]=(j*i+(1-j)*l);
		temp[1]=j;;//Forward=1, Reverse Complement=0
		temp[2]=(j*l+(1-j)*i);
		temp[3]=1-j;;//Forward=1, Reverse Complement=0
		temp[4]=unitigs[i][j*unitigs_info[i][0]][4+k*5];
						  
		unitigs_con[i][j].push_back(temp);

		
		temp[0]=(j*l+(1-j)*i);
		temp[1]=1-j;
		temp[2]=(j*i+(1-j)*l);
		temp[3]=j;
		temp[4]=unitigs[i][j*unitigs_info[i][0]][4+k*5];

		unitigs_con[l][m].push_back(temp);
	      }
	    }
	  }
	}
      }
    }
  }
}

void check_self_pairing(vector<vector<vector<int> > > &unis,int mate_table[num_of_reads+1],int &num_of_unitigs,vector<vector<int> > &unitigs_info){

  int i,j,k;
  int mate_pair_distance;

   for (i=0;i<num_of_unitigs;i++){
    for (j=0;j<unitigs_info[i][0]-1;j++){ 
      for (k=j+1;k<unitigs_info[i][0];k++){
	if( ((unis[i][j][0]==unis[i][k][0]+1) && (unis[i][k][0]%2==0))  ||  ((unis[i][j][0]==unis[i][k][0]-1) && (unis[i][k][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis[i][j][1]==1 && unis[i][k][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	    
	    mate_pair_distance=unis[i][k][2]-unis[i][j][2];

	    if ((l_bd_mp <=mate_pair_distance ) && (u_bd_mp >=mate_pair_distance )){//two reads are distanced between 2400-3600
		mate_table[unis[i][j][0]]=1;mate_table[unis[i][k][0]]=1;//found mate-pair within a unitig
		mate_table[num_of_reads]+=2;
	    }
	  }
	}
      }
    }
   }

#if 0
   for (i=0;i<num_of_reads+1;i++){
     cout<<i<<"  "<<mate_table[i]<<"\n";
   }
#endif
   
}

/*
void mate_pair_check(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info,int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con, vector<vector<int> > &contig_unis_list, int mate_table[num_of_reads+1]){

  int i,j,k,l,m,n;
  int FB,next_uni;
  int mate_pair_count, fatal_no_mate_count,mate_pair_distance;

  
  i=contig_unis_list[(int)(contig_unis_list.size())-1][0];
  for (j=0;j<2;j++){ 
    for (k=0;k<unitigs_con_count[i][j+1];k++){
	
      FB=unitigs_con[i][j][k][(1-j) + j*3];//1= Forward or 0= Reverse complement
      next_uni=unitigs_con[i][j][k][j*2];// the index for compared unitig
      mate_pair_count=0;
      fatal_no_mate_count=0;
	
      if(FB==1){//connected forward
	for (l=0;l<unitigs_info[i][0];l++){
	  for (m=0;m<unitigs_info[next_uni][0];m++){
	    if( ((unis[i][l][0]==unis[next_uni][m][0]+1) && (unis[next_uni][m][0]%2==0))  || ((unis[i][l][0]==unis[next_uni][m][0]-1) && (unis[next_uni][m][0]%2==1))){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	      if(unis[i][l][1]==j && unis[next_uni][m][1]==(1-j) ){// two reads are facing each other 5'-3' 3'-5' way

		// cout<<unis[next_uni][m][0]<<"  " <<unis[i][l][0]<<"  ";


		  
		mate_pair_distance=(1-j)*(unitigs_info[next_uni][1]-unis[next_uni][m][2]+unis[i][l][2])+j*(unitigs_info[i][1]-unis[i][l][2]+unis[next_uni][m][2]);

		//		  cout<<mate_pair_distance<<"  ";


		if ((l_bd_mp <=mate_pair_distance ) && (u_bd_mp >=mate_pair_distance )){//two reads are distanced between 2400-3600
		  mate_pair_count++;
		}
	      }
	    }
	  }
	}
      }
      else {//connected reverse_complement
	for (l=0;l<unitigs_info[i][0];l++){
	  for (m=0;m<unitigs_info[next_uni][0];m++){
	    if( ((unis[i][l][0]==unis_RC[next_uni][m][0]+1) && (unis_RC[next_uni][m][0]%2==0))  || ((unis[i][l][0]==unis_RC[next_uni][m][0]-1) && (unis_RC[next_uni][m][0]%2==1))){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	      if(unis[i][l][1]==j && unis_RC[next_uni][m][1]==(1-j) ){// two reads are facing each other 5'-3' 3'-5' way
		 
		mate_pair_distance=(1-j)*(unitigs_info[next_uni][1]-unis_RC[next_uni][m][2]+unis[i][l][2])+j*(unitigs_info[i][1]-unis[i][l][2]+unis_RC[next_uni][m][2]);
		if ((l_bd_mp <=mate_pair_distance ) && (u_bd_mp >=mate_pair_distance )){//two reads are distanced between 2400-3600
		  mate_pair_count++;
		}
	      }
	    }
	  }
	}
      }
	
      //	cout<<i<<"  "<<j<<"  "<<next_uni<<"  "<<unitigs_con[i][j][k][(1-j) + j*3]<<"  "<<mate_pair_count<<"\n";
	
    }
  }
}
*/

int mate_pair_check(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC, vector<vector<int> > &unitigs_info, int unitigs_con_count[][3], int mate_table[num_of_reads+1], int uni1, int uni1_FR, int uni2, int uni2_FR, int distance){

  int i,j,mate_pair_distance;
  int mate_pair_count=0;
  
  for (i=0;i<unitigs_info[uni1][0];i++){
    for (j=0;j<unitigs_info[uni2][0];j++){
      if (uni1_FR==1 && uni2_FR==1){
	if( ((unis[uni1][i][0]==unis[uni2][j][0]+1) && (unis[uni2][j][0]%2==0)) || ((unis[uni1][i][0]==unis[uni2][j][0]-1) && (unis[uni2][j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis[uni1][i][1]==1 && unis[uni2][j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  
	    mate_pair_distance=distance+unitigs_info[uni1][1]-unis[uni1][i][2]+unis[uni2][j][2];
	    if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) ){//two reads are distanced between 2400-3600
	      mate_table[num_of_reads]+=2;
	      mate_table[unis[uni1][i][0]]+=1;mate_table[unis[uni2][j][0]]+=1;
	      mate_pair_count+=1;
	    }	  
	  }	  
	}
      }
      else if (uni1_FR==1 && uni2_FR==0){
	if( ((unis[uni1][i][0]==unis_RC[uni2][j][0]+1) && (unis_RC[uni2][j][0]%2==0)) || ((unis[uni1][i][0]==unis_RC[uni2][j][0]-1) && (unis_RC[uni2][j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis[uni1][i][1]==1 && unis_RC[uni2][j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  
	    mate_pair_distance=distance+unitigs_info[uni1][1]-unis[uni1][i][2]+unis_RC[uni2][j][2];
	    if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) ){//two reads are distanced between 2400-3600
	      mate_table[num_of_reads]+=2;
	      mate_table[unis[uni1][i][0]]+=1;mate_table[unis_RC[uni2][j][0]]+=1;
	      mate_pair_count+=1;
	    }	  
	  }	  
	}
      }   
      else if (uni1_FR==0 && uni2_FR==1){
	if( ((unis_RC[uni1][i][0]==unis[uni2][j][0]+1) && (unis[uni2][j][0]%2==0)) || ((unis_RC[uni1][i][0]==unis[uni2][j][0]-1) && (unis[uni2][j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis_RC[uni1][i][1]==1 && unis[uni2][j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  
	    mate_pair_distance=distance+unitigs_info[uni1][1]-unis_RC[uni1][i][2]+unis[uni2][j][2];
	    if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) ){//two reads are distanced between 2400-3600
	      mate_table[num_of_reads]+=2;
	      mate_table[unis_RC[uni1][i][0]]+=1;mate_table[unis[uni2][j][0]]+=1;
	      mate_pair_count+=1;
	    }	  
	  }	  
	}
      }
      else{ //uni1_FR==0 && uni2_FR==0
	if( ((unis_RC[uni1][i][0]==unis_RC[uni2][j][0]+1) && (unis_RC[uni2][j][0]%2==0)) || ((unis_RC[uni1][i][0]==unis_RC[uni2][j][0]-1) && (unis_RC[uni2][j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis_RC[uni1][i][1]==1 && unis_RC[uni2][j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  
	    mate_pair_distance=distance+unitigs_info[uni1][1]-unis_RC[uni1][i][2]+unis_RC[uni2][j][2];
	    if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) ){//two reads are distanced between 2400-3600
	      mate_table[num_of_reads]+=2;
	      mate_table[unis_RC[uni1][i][0]]+=1;mate_table[unis_RC[uni2][j][0]]+=1;
	      mate_pair_count+=1;
	    }	  
	  }	  
	}
      }
    }
  }
  
  return mate_pair_count;
}

int iterate_for_finding_a_contig(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info,int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con, vector<vector<int> > &contig_unis_list, int mate_table[num_of_reads+1]){

  int contig_size=(int)(contig_unis_list.size());
  int i,j,k;
  int uni1,uni1_FR,uni2,uni2_FR;
  int last_uni, last_uni_FR;
  
  int distance, distance0;
  int mate_pair_count;
  
  vector<vector<int> > next_uni_table; //table of next possible unitigs

  vector<int> temp(4);
  
  last_uni=contig_unis_list[contig_size-1][0];
  last_uni_FR=contig_unis_list[contig_size-1][1];

  for (k=0;k<unitigs_con_count[last_uni][last_uni_FR+1];k++){

    uni2=unitigs_con[last_uni][last_uni_FR][k][last_uni_FR*2];
    uni2_FR=last_uni_FR*unitigs_con[last_uni][last_uni_FR][k][last_uni_FR*2+1]+(1-last_uni_FR)*(1-unitigs_con[last_uni][last_uni_FR][k][last_uni_FR*2+1]);
    mate_pair_count=0;
    
    for(i=0;i<contig_size;i++){

      j=contig_size-i-1;
    
      uni1=contig_unis_list[j][0];
      uni1_FR=contig_unis_list[j][1];
  
      if (i==0)
	distance0=0;
      else
	distance0=distance0+unitigs_info[contig_unis_list[j+1][0]][1]+contig_unis_list[j+1][2];
    
      distance=distance0+unitigs_con[last_uni][last_uni_FR][k][4];
  
      if (distance<=u_bd_mp && distance>=l_bd_mp ){
	mate_pair_count+= mate_pair_check(unis, unis_RC, unitigs_info, unitigs_con_count, mate_table, uni1, uni1_FR, uni2, uni2_FR, distance);
      }
    }

    temp[0]=uni2;
    temp[1]=uni2_FR;
    temp[2]=-read_len+unitigs_con[last_uni][last_uni_FR][k][4];//distance between contig and the next unitig(negative number of course)
    temp[3]=mate_pair_count;

    cout<<temp[3]<<" ";
    
    next_uni_table.push_back(temp);

  }
  // iterate_for_finding_a_contig(unis, unis_RC,num_of_unitigs,unitigs_info,unitigs_con_count, unitigs_con,contig_unis_list, mate_table);
}

void really_find_a_contig(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info,int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con, char *contig){
 
  int i,start_unitig=0;
  
  for (i=0;i<num_of_unitigs;i++){//if an end of a unitig does not have a connecting edges, then the unitig can be a boundary of a contig
    if(unitigs_con_count[i][1]==0 || unitigs_con_count[i][2]==0){
      if(unitigs_con_count[i][0]==0){
	cout<<"There is a unitig without any edges to outside";
      }
      else{
	start_unitig=i;
	break;
      }
    }
  }

  vector<vector<int> > contig_unis_list; //unitig #, F/RC, starting pt
  vector<int> a_unis(3);
  a_unis[0]=start_unitig;
  a_unis[1]=1;
  a_unis[2]=0;
  contig_unis_list.push_back(a_unis);//insert the first unitig

  int mate_table[num_of_reads+1]={0};//count mate_pair numbers
  check_self_pairing(unis,mate_table,num_of_unitigs,unitigs_info);
  
  iterate_for_finding_a_contig(unis, unis_RC, num_of_unitigs, unitigs_info, unitigs_con_count, unitigs_con,contig_unis_list,mate_table);

}

void  find_a_contig(vector<vector<vector<int> > > &unitigs,vector<vector<int> > &unitigs_info, char *contig){

  int num_of_unitigs=unitigs.size();//number of unitigs
  int unitigs_con_count[num_of_unitigs][3];//number of connected unitigs count
  int num_of_connections=count_the_num_of_connections(unitigs, num_of_unitigs, unitigs_info, unitigs_con_count);//number of total connections
  
  vector<vector<vector<int> > > unis(num_of_unitigs);
  vector<vector<vector<int> > > unis_RC(num_of_unitigs);
  setup_unis(unitigs, unis, unis_RC,num_of_unitigs, unitigs_info );
  
  int i,j,k;
 
  vector<vector<vector<vector<int> > > > unitigs_con(num_of_unitigs,vector<vector<vector<int> > >(2)); //record connected unitigs [prior unitig F/B next unitig F/B distance]
  connected_unitigs(unitigs, num_of_unitigs,unitigs_info, unitigs_con_count,unitigs_con);

  really_find_a_contig(unis,unis_RC,num_of_unitigs,unitigs_info,unitigs_con_count,unitigs_con, contig);





#if 0
  cout<<"unitigs_info"<<"\n";
  for (i=0;i<num_of_unitigs;i++)
    {
      for(j=0;j<2;j++){
	cout<<unitigs_info[i][j]<<"   ";
      }
      cout<<"\n";
    }

   cout<<"unitigs_con_count"<<"\n";
  for (i=0;i<num_of_unitigs;i++){
      for(j=0;j<3;j++){
	cout<<unitigs_con_count[i][j]<<"   ";
      }
      cout<<"\n";
    }
    
  cout<<"unitigs_con"<<"\n";
  for (i=0;i<num_of_unitigs;i++){
      for(j=0;j<2;j++){
	cout<<i<<" "<<j<<" ";
	for(k=0;k<unitigs_con[i][j].size();k++){

	  for (int l=0;l<unitigs_con[i][j][k].size();l++)
	   cout<<unitigs_con[i][j][k][l]<<"   ";

	}
	 cout<<"\n";
  
      }
     
    }
  
   for (i=0;i<unitigs_info[2][0];i++){
     cout<< unis[2][i][0]<<" "<<unis[2][i][1]<<" "<<unis[2][i][2]<<"     "<< unis_RC[2][i][0]<<" "<< unis_RC[2][i][1]<<" "<< unis_RC[2][i][2]<<"\n";
   }
#endif


}

int main(){

  int i,j,k;
  
  //read from fasta file
  char list_of_reads[num_of_reads][read_len + 1];
  for (i=0;i<num_of_reads;i++){
    list_of_reads[i][0] = '\0';
  }

  read_from_fasta(list_of_reads);

  //define reverse complement reads
  char  list_of_reads_RC[num_of_reads][read_len + 1];
  for (i=0;i<num_of_reads;i++){
    list_of_reads_RC[i][0] = '\0';
  }
  
  reversecomplement(list_of_reads, list_of_reads_RC);

  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][6];//read1,read2,F/R(1/0),olap_length,first_read_location F/B(1/0) , deleted(1/0)

  vector<vector<int> > list_of_exact_olaps;//save exact overlaps
  int num_of_exact_olaps=0;
  int num_of_edges_to_delete=0;
  int location[num_of_reads+1]={0};//starting point of the read 1 in the list_of_olaps

  int num_of_olaps=read_from_olaps(list_of_olaps,location,list_of_exact_olaps,num_of_exact_olaps,num_of_edges_to_delete);

  record_edge_to_delete(list_of_olaps,num_of_olaps,location, num_of_edges_to_delete);

  int num_of_edges_for_unitigs=num_of_olaps-num_of_edges_to_delete;

  vector<vector<vector<int> > > edges_for_nodes(num_of_reads, vector<vector<int> >(2));//for each node, save outgoing edges, incoming edges
  int edges_for_nodes_index[num_of_reads][4]={{0}}; //save #of total edges, # of outgoing edges,  # of incoming edges, #used or not
  set_up_viable_edges(location, list_of_olaps, edges_for_nodes,edges_for_nodes_index, list_of_exact_olaps, num_of_exact_olaps);
  vector<vector<vector<int> > > edges_for_nodes_RC(num_of_reads, vector<vector<int> >(2));//reverse complement list
  set_up_edges_RC(edges_for_nodes, edges_for_nodes_RC);//setup RC
 
  vector<vector<vector<int> > > unitigs;// contains all unitigs: unitig_number,edges 
  vector<vector<int> > unitigs_info; //contains # of reads, total lengths of unitigs
  find_unitigs(unitigs,unitigs_info,edges_for_nodes,edges_for_nodes_RC, edges_for_nodes_index,list_of_exact_olaps, num_of_exact_olaps);

  char contig[num_of_reads*read_len];
  find_a_contig(unitigs,unitigs_info,contig);

  //print_a_contig(contig);
  
  return 0;
}




 

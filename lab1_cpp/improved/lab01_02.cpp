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


//
//HW2 codes for finding unitigs start here
//

class read_2{
public:
  int num;//=read number 0-299
  bool direction;//1=forward 0=reverse complement
  read_2 rc();
};

read_2 read_2::rc(){
  read_2 temp;
  temp.num=num;
  temp.direction=1-direction;
  return temp;
}

class edge_b {
public:
  read_2 from_read, to_read;//from_read to to_read
  edge_b rc();
};

edge_b edge_b::rc(){
  edge_b temp;
  temp.from_read=to_read.rc();
  temp.to_read=from_read.rc();
  return temp;
}

class edge {
  edge_b olap_f, olap_rc;
  int offset;//how read 1 is ahead of read 2
  bool deleted;//deleted=1 if deleted edge
  void left_to_right_align_and_set_rc(); //change the directions of edges to left to right and set up reverse_complement edge
};

void edge::left_to_right_align_and_set_rc(){
  read_2 temp;
  if (offset<0){//left_to_right_align
    offset=-offset;
    temp=olap_f.from_read;
    olap_f.from_read=olap_f.to_read;
    olap_f.to_read=temp;
  }
  olap_rc=olap_f.rc();
}

/*




list_of_olaps[num_of_reads*(num_of_reads-1)/2][6];//read1,read2,F/R(1/0),olap_length,first_read_location F/B(1/0) , deleted(1/0)





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

*/

//
//HW2 codes for finding unitigs end here
//


int main(){

  //HW2 finds unitigs

  //HW2 ends

  /*
  int list_of_olaps[num_of_reads*(num_of_reads-1)/2][6];//read1,read2,F/R(1/0),olap_length,first_read_location F/B(1/0) , deleted(1/0)

  vector<vector<int> > list_of_exact_olaps;//save exact overlaps
  int num_of_exact_olaps=0;
  int num_of_edges_to_delete=0;
  int location[num_of_reads+1]={0};//starting point of the read 1 in the list_of_olaps

  int num_of_olaps=read_from_olaps(list_of_olaps,location,list_of_exact_olaps,num_of_exact_olaps,num_of_edges_to_delete);


  record_edge_to_delete(list_of_olaps,num_of_olaps,location, num_of_edges_to_delete);
  */ 
 return 0;
}



  /*
  
  int num_of_edges_for_unitigs=num_of_olaps-num_of_edges_to_delete;

  vector<vector<vector<int> > > edges_for_nodes(num_of_reads, vector<vector<int> >(2));//for each node, save outgoing edges, incoming edges
  int edges_for_nodes_index[num_of_reads][4]={{0}}; //save #of total edges, # of outgoing edges,  # of incoming edges, #used or not
  set_up_viable_edges(location, list_of_olaps, edges_for_nodes,edges_for_nodes_index, list_of_exact_olaps, num_of_exact_olaps);
  vector<vector<vector<int> > > edges_for_nodes_RC(num_of_reads, vector<vector<int> >(2));//reverse complement list
  set_up_edges_RC(edges_for_nodes, edges_for_nodes_RC);//setup RC
 
  vector<vector<vector<int> > > unitigs;// contains all unitigs: unitig_number,edges
  vector<vector<int> > unitigs_info; //contains # of reads, total lengths of unitigs 
  find_unitigs(unitigs,unitigs_info,edges_for_nodes,edges_for_nodes_RC, edges_for_nodes_index,list_of_exact_olaps, num_of_exact_olaps);

  print_unis(unitigs,unitigs_info);
  */
 



/*

void record_edge_to_delete(int list_of_olaps[][6],int num_of_olaps, int location[num_of_reads+1],int &num_of_edges_to_delete){

//i the number of read1, j the location of read1, read2 pair, k the location of read1, read2 pair, l, the number of
		read2, read2 pair

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

void print_unis( vector<vector<vector<int> > > &unitigs,vector<vector<int> > &unitigs_info){

  int i,j,k,ind,sum;
  FILE * pFile;

#if SAMPLE
  pFile = fopen ("sample.unis","w");
#else
  pFile = fopen ("lab01.unis","w");
#endif
 	
	
  for (i=0;i<unitigs.size();i++){
 
    fprintf (pFile, "UNI  %02d %*d %*d\n", i+1, 5,unitigs_info[i][0],6,unitigs_info[i][1]);//print the title

    fprintf (pFile, "  %03d  ",unitigs[i][1][0]+1);//print the first element
    ind=unitigs[i][1][1];
    if (ind==1)
      fprintf (pFile, "F  ");
    else
      fprintf (pFile, "R  ");
    fprintf (pFile, "%*d\n",4,0);


    for (j=1;j<unitigs_info[i][0];j++){//print the rest
      
      fprintf (pFile, "  %03d  ",unitigs[i][j][2]+1);
      ind=unitigs[i][j][3];
      if (ind==1)
	fprintf (pFile, "F  ");
      else
	fprintf (pFile, "R  ");
      fprintf (pFile, "%*d\n",4,unitigs[i][j][4]);
    }
  }
}



void print_unis_raw ( vector<vector<vector<int> > > &unitigs){

  int i,j,k;
  FILE * pFile;

#if SAMPLE
  pFile = fopen ("sample.unis_raw","w");
#else
  pFile = fopen ("lab01.unis_raw","w");
#endif
 	

  for (i=0;i<unitigs.size();i++){
    for (j=0;j<unitigs[i].size();j++){
      for (k=0;k<unitigs[i][j].size();k++){
	fprintf (pFile, "%d  ",unitigs[i][j][k]);
      }
      fprintf (pFile, "\n");
    }
    fprintf (pFile, "\n\n");
  }
}

int main(){

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

  print_unis(unitigs,unitigs_info);
  
  return 0;
}
*/

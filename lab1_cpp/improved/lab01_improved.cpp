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
//HW1 codes for finding overlapss start here
//

class read_raw{
public:
  char read[read_len+1]; //read
  char read_rc[read_len+1]; //reverse complement of the read
  read_raw();
  void reverse_complement();
};

read_raw::read_raw () {
  read[0]='\0';
}

void read_raw::reverse_complement(){
  for (int i=0;i<read_len;i++){
    if (read[read_len-i-1]=='a')
      read_rc[i]='t';
    else if (read[read_len-i-1]=='t')
      read_rc[i]='a'; 
    else if (read[read_len-i-1]=='c')
      read_rc[i]='g';
    else 
      read_rc[i]='c';
  }
  read_rc[read_len]='\0';
}

class olap {
public:
  int f_read,t_read; //from_read and to_read
  bool ori_t; //  ori_t=1/0 if t_read is forward/reverse complement
  int offset;
};

class olaps {
public:
  int size;
  vector<olap> list;
};

class read_end_piece{
public:
  char end[41];
  int failure[40];
  bool location;
  read_end_piece();
};

read_end_piece::read_end_piece () {
  end[40]='\0';
}

class read_ends{
public:
  read_end_piece right, left;
  read_ends();
};

read_ends::read_ends () {
  right.location=0;
  left.location=1;
}

void read_from_fasta(read_raw list_of_reads[num_of_reads]){

  char str[read_len];
  int i=0;

#if SAMPLE
  ifstream infile("sample.fasta");
#else
  ifstream infile("lab01.fasta");
#endif

  if (infile.is_open()){
    while (infile >> str){
      if (str[0]=='>'){
	i=atoi(str+1)-1;
      }
      else{
	strcat(list_of_reads[i].read,str);
      }
    }
  }
  infile.close();
}

void reverse_complement(read_raw list_of_reads[num_of_reads]){
  for (int i=0;i<num_of_reads;i++){
    list_of_reads[i].reverse_complement();
  }
}

void KMP_table(read_end_piece &piece){

  int cnd;
  piece.failure[0] = -1;

  for (int pos=1;pos < 40;pos++){
    cnd=piece.failure[pos-1];

    while((piece.end[pos]!=piece.end[cnd+1]) && (cnd>=0)){
      cnd=piece.failure[cnd];
    }

    if (piece.end[pos]==piece.end[cnd+1]){
      piece.failure[pos]=cnd+1;
    }
    else {
      piece.failure[pos]=-1;
    }
  }
}

bool check_rest(bool LEFT, int match_start, int match_end, char *first_read, char *second_read){
  int i;
  if (LEFT){//left
    for(i=match_end+1;i<read_len;i++){
      if (first_read[i-match_start]!=second_read[i])
	return 0;
    }
  }
  else {//right
    for (i=0;i<match_start;i++){
      if (first_read[read_len-match_end+i-1]!=second_read[i])
	return 0;
    }
  }
  return 1;
}

bool KMP_search(char *first_read, read_end_piece &first_read_piece, char *second_read, olap &Olap){

  int i=0,j,matched_len=0, match_start,match_end;

  while (i < read_len){
    if (first_read_piece.end[matched_len] == second_read[i]){
      i++;
      matched_len++;
    }
    else if (matched_len==0)
      i++;
    else
      matched_len=first_read_piece.failure[matched_len-1]+1;
    if (matched_len==40){
      match_start=i-40;
      match_end=i-1;
      if (check_rest(first_read_piece.location, match_start,match_end,first_read,second_read)){
	Olap.offset=-first_read_piece.location*match_start+(1-first_read_piece.location)*(read_len-match_end-1);
	return 1;
      }
    }
  }
  return 0;
}

bool find_overlaps_of_two_reads(read_raw &first, read_raw &second,olap &Olap, read_ends &Read_ends){

  bool found=0;
  
  if (KMP_search(first.read, Read_ends.right, second.read, Olap)){
    Olap.ori_t=1;
    found=1;
  }
  else if (KMP_search(first.read, Read_ends.right, second.read_rc, Olap)){
    Olap.ori_t=0;
    found=1;
  }
  else if (KMP_search(first.read, Read_ends.left, second.read, Olap)){
    Olap.ori_t=1;
    found=1;
  }
  else if (KMP_search(first.read, Read_ends.left, second.read_rc, Olap)){
    Olap.ori_t=0;
    found=1;
  }
  
  return found;
}

void find_olaps(read_raw list_of_reads[num_of_reads], olaps &list_of_olaps){

  int first_read, second_read;
  olap Olap;
  read_ends Read_ends;
  
  list_of_olaps.size=0;
  
  for (first_read=0;first_read<num_of_reads-1;first_read++){

    strncpy(Read_ends.right.end, list_of_reads[first_read].read+read_len-40,40);
    KMP_table(Read_ends.right);

    strncpy(Read_ends.left.end, list_of_reads[first_read].read,40);
    KMP_table(Read_ends.left);    	

    for (second_read=first_read+1;second_read<num_of_reads;second_read++){
      Olap.f_read=first_read;
      Olap.t_read=second_read;
      if (find_overlaps_of_two_reads(list_of_reads[first_read],list_of_reads[second_read],Olap,Read_ends)){
	list_of_olaps.size++;
	list_of_olaps.list.push_back(Olap);
      }
    }
  }
}
 
void print_olaps(olaps &list_of_olaps){

  FILE * pFile;

#if SAMPLE
  pFile = fopen ("sample.olaps","w");
#else
  pFile = fopen ("lab01.olaps","w");
#endif

  for (int i=0;i<list_of_olaps.size;i++){
    fprintf (pFile, " %03d  ",list_of_olaps.list[i].f_read+1);
    fprintf (pFile, "%03d  ",list_of_olaps.list[i].t_read+1);
    if (list_of_olaps.list[i].ori_t)
      fprintf (pFile, "F  ");
    else
      fprintf (pFile, "R  ");
    fprintf (pFile, "%*d\n",4,list_of_olaps.list[i].offset);
  }
  
  fclose (pFile);
}

void find_and_print_olaps(read_raw list_of_reads[num_of_reads], olaps &list_of_olaps){

  read_from_fasta(list_of_reads); 
  reverse_complement(list_of_reads);
  find_olaps(list_of_reads,list_of_olaps);
  print_olaps(list_of_olaps);
}

//
//HW1 codes for finding overlaps end here
//



//
//HW2 codes for finding unitigs start here
//
class node {
public:
  int num;//=read number 0-299
  bool ori;//1=forward 0=reverse complement
};

void rc_node(node &a){//reverse complement a read
  a.ori=1-a.ori;
}

class edges_node {
public:
  vector<node> f_node;//nodes from which edges come to this node
  vector<int> f_offset;//how much from node is ahead of to node
  
  vector<node> t_node;//nodes to which edges go from this node
  vector<int> t_offset;//how much current node is ahead of to node

  int total_edge_ct;// #of total edges connected to the node
  int f_edge_ct;// # of edges coming to the node
  int t_edge_ct;//# of edges going out from the node
  bool used;//used=1 if the node is used
  edges_node();
};

edges_node::edges_node(){
  total_edge_ct=0;
  f_edge_ct=0;
  t_edge_ct=0;
  used=0;
}

class olap_2: public olap {
public:
  bool r_arrow; //direction of arrow is right
  bool deleted; //deleted=1
};

class olaps_2 {
public:  
  int size;//total size
  int deleted_size;//deleted overlap size
  int exact_size;//exact overlap size
  vector<olap_2> list;//overlap list
  int f_read_loc[num_of_reads+1];//starting point of the f_read in the list_of_olaps
  vector<olap_2> exact;//save the exact overlaps in the list_of_olaps

  edges_node nodes[num_of_reads];
  olaps_2();
};

olaps_2::olaps_2(){
  size=0; deleted_size=0; exact_size=0;
  fill(f_read_loc, f_read_loc+num_of_reads+1,0);
}

class unitigs{
public:
  int size;
  //  int num_of_edges_for_unitigs=num_of_olaps-num_of_edges_to_delete;
  //vector<vector<vector<int> > > unitigs;// contains all unitigs: unitig_number,edges
  //  vector<vector<int> > unitigs_info; //contains # of reads, total lengths of unitigs
};













void read_from_olaps(olaps_2 &l_olaps){

  char read_1_c[5],read_2_c[5], read_2_FR_c[5], offset_c[5];//read 1, read 2, Front/RC of read 2, offset
  int i,j,offset,read_1,num_of_olaps=0;
  olap_2 temp_olap;

#if SAMPLE
  ifstream infile("sample.olaps");
#else
  ifstream infile("lab01.olaps");
#endif
  
  if (infile.is_open()) {
    while (infile >> read_1_c>> read_2_c>> read_2_FR_c >> offset_c){

      read_1=atoi(read_1_c);

      temp_olap.f_read=read_1-1;
      temp_olap.t_read=atoi(read_2_c)-1;

      if (read_2_FR_c[0]=='F')
	temp_olap.ori_t=1;
      else
	temp_olap.ori_t=0;

      offset=atoi(offset_c);

      if( offset==0){ // two reads are identical
	l_olaps.exact_size++;
	l_olaps.exact.push_back(temp_olap);
      }
	
      else{//two reads are not identical
       	
	temp_olap.offset=offset;
	if ( offset >0)
	  temp_olap.r_arrow=1;
	else 
	  temp_olap.r_arrow=0;	
	l_olaps.f_read_loc[read_1]=num_of_olaps;
	temp_olap.deleted=0;
	l_olaps.list.push_back(temp_olap);
	num_of_olaps++;
      }
    }
    infile.close();
  }
  
  for (i=1;i<num_of_reads;i++){//fill the f_read location table
    if (l_olaps.f_read_loc[i]==0)
      l_olaps.f_read_loc[i]=l_olaps.f_read_loc[i-1];
    else l_olaps.f_read_loc[i]+=1;
  }
  l_olaps.f_read_loc[num_of_reads]=num_of_olaps;

  for (i=0;i<l_olaps.exact_size;i++){//delete olaps connecting exact olaps
    for (j=0;j<num_of_olaps;j++){
      if((l_olaps.exact[i].t_read==l_olaps.list[j].f_read || l_olaps.exact[i].t_read==l_olaps.list[j].t_read) && l_olaps.list[j].deleted==0){
	l_olaps.list[j].deleted=1;
	l_olaps.deleted_size++;
      }
    }
  }

  l_olaps.size=num_of_olaps;
}

void record_edge_to_delete(olaps_2 &l_olaps){

//i the number of read1
//j the location of read1, read2_1 pair
//k the location of read1, read2_2 pair
//l the number of read2_1, read2_2 pair

  int i,j,k,l;
  bool deleted=1;
  int arrow0,arrow1,arrow2;
   
  for (i=0;i<num_of_reads-2;i++){
    for (j=l_olaps.f_read_loc[i];j<l_olaps.f_read_loc[i+1]-1;j++){
      for (k=j+1;k<l_olaps.f_read_loc[i+1];k++){
	for (l=l_olaps.f_read_loc[l_olaps.list[j].t_read];l<l_olaps.f_read_loc[l_olaps.list[j].t_read+1];l++){
      	  if(l_olaps.list[k].t_read==l_olaps.list[l].t_read){
    	    if((l_olaps.list[j].ori_t+l_olaps.list[k].ori_t+l_olaps.list[l].ori_t)%2==1){
	      arrow0=l_olaps.list[j].r_arrow;
	      arrow1=l_olaps.list[k].r_arrow;
	      arrow2=l_olaps.list[l].r_arrow;
	      if (l_olaps.list[j].ori_t==0)
		arrow2=1-arrow2;
	      if (arrow0==arrow1 && arrow1==arrow2 && l_olaps.list[k].deleted==0){
		l_olaps.list[k].deleted=deleted;
		l_olaps.deleted_size++;
	      }
	      else if (arrow0==arrow1 && arrow0!=arrow2 && l_olaps.list[j].deleted==0){
		l_olaps.list[j].deleted=deleted;
		l_olaps.deleted_size++;
	      }
	      else if (arrow0!=arrow1 && arrow1==arrow2 && l_olaps.list[l].deleted==0){
    		l_olaps.list[l].deleted=deleted;
		l_olaps.deleted_size++;
      	      }
	    }
	  }
	  else if (l_olaps.list[l].t_read>l_olaps.list[k].t_read)
	    break;
	}
      }
    }
  }
}


void set_up_viable_edges(olaps_2 &l_olaps){

  int i,j;
  node temp_node;
  
  for (i=0;i<num_of_reads;i++){
    for (j=l_olaps.f_read_loc[i];j<l_olaps.f_read_loc[i+1];j++){
      if (l_olaps.list[j].deleted==0){ //not deleted
	if (l_olaps.list[j].r_arrow){ //location of first read front

	  temp_node.num=l_olaps.list[j].t_read;
	  temp_node.ori=l_olaps.list[j].ori_t;
	  l_olaps.nodes[i].t_node.push_back(temp_node); 

	  l_olaps.nodes[i].t_offset.push_back(l_olaps.list[j].offset);

	  l_olaps.nodes[i].total_edge_ct+=1;
	  l_olaps.nodes[i].t_edge_ct+=1;

	  if (l_olaps.list[j].ori_t){//second read is front oriented

	    temp_node.num=i;
	    temp_node.ori=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].f_node.push_back(temp_node); 
	    
	    l_olaps.nodes[l_olaps.list[j].t_read].f_offset.push_back(l_olaps.list[j].offset);

	    l_olaps.nodes[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].f_edge_ct+=1;	    
	  }
	  else{//second read is reverse complement oriented
	 
	    temp_node.num=i;
	    temp_node.ori=0;
	    l_olaps.nodes[l_olaps.list[j].t_read].t_node.push_back(temp_node); 
    	
	    l_olaps.nodes[l_olaps.list[j].t_read].t_offset.push_back(l_olaps.list[j].offset);

	    l_olaps.nodes[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].t_edge_ct+=1;	    
	  }
	}
	else { //location of first read back
	
	  temp_node.num=l_olaps.list[j].t_read;
	  temp_node.ori=l_olaps.list[j].ori_t;
	  l_olaps.nodes[i].f_node.push_back(temp_node); 
	
	  l_olaps.nodes[i].f_offset.push_back(-l_olaps.list[j].offset);

	  l_olaps.nodes[i].total_edge_ct+=1;
	  l_olaps.nodes[i].f_edge_ct+=1;

	  if (l_olaps.list[j].ori_t){//second read is front oriented
 
	    temp_node.num=i;
	    temp_node.ori=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].t_node.push_back(temp_node); 

	    l_olaps.nodes[l_olaps.list[j].t_read].t_offset.push_back(-l_olaps.list[j].offset);

	    l_olaps.nodes[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].t_edge_ct+=1;	    
	  }
	    
	  
	  else{//second read is reverse complement oriented
	  	   
	    temp_node.num=i;
	    temp_node.ori=0;
	    l_olaps.nodes[l_olaps.list[j].t_read].f_node.push_back(temp_node); 
    	
	    l_olaps.nodes[l_olaps.list[j].t_read].f_offset.push_back(-l_olaps.list[j].offset);

	    l_olaps.nodes[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.nodes[l_olaps.list[j].t_read].f_edge_ct+=1;	    
	  }
	}
      }
    }
  }

  for (i=0;i<l_olaps.exact_size;i++){
    l_olaps.nodes[l_olaps.exact[i].t_read].used=1;
  }
}


void find_a_unitig(int &starting_point, edges_node e_for_nodes[num_of_reads], unitigs &unis){

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


void find_unitigs(edges_node e_for_nodes[num_of_reads], unitigs &unis){
  
  int read_exact_match_count=0;
  
  for(int starting_point=0;starting_point<num_of_reads;starting_point++){
    if (!e_for_nodes[starting_point].used){//node not used
      e_for_nodes[starting_point].used=1;
      find_a_unitig(starting_point, e_for_nodes, unis);
    }
  }
}   

void find_unis(unitigs &unis){

  olaps_2 l_olaps; 
  edges_node e_for_nodes[num_of_reads];
  read_from_olaps(l_olaps);
  record_edge_to_delete(l_olaps);
  set_up_viable_edges(l_olaps,e_for_nodes);
  // set_up_edges_RC(edges_for_nodes);
  find_unitigs(e_for_nodes, unis);
 

#if 1
  FILE * pFile;
  pFile = fopen ("lab01.temp","w");

  
  for (int i=0;i<num_of_reads;i++){
    fprintf (pFile, "%d %d %d %d ",i,e_for_nodes[i].total_edge_ct,e_for_nodes[i].f_edge_ct,e_for_nodes[i].t_edge_ct);
    for (int j=0;j<e_for_nodes[i].f_edge_ct;j++)
      fprintf (pFile, "%d %d %d 1 %d ",  e_for_nodes[i].f_node[j].num,e_for_nodes[i].f_node[j].ori ,i,e_for_nodes[i].f_offset[j]);
    fprintf (pFile, "\n");
    fprintf (pFile, "%d %d %d %d ",i,e_for_nodes[i].total_edge_ct,e_for_nodes[i].f_edge_ct,e_for_nodes[i].t_edge_ct);

    for (int j=0;j<e_for_nodes[i].t_edge_ct;j++)
      fprintf (pFile, "%d 1 %d %d %d ",i, e_for_nodes[i].t_node[j].num,e_for_nodes[i].t_node[j].ori, e_for_nodes[i].t_offset[j]);
    fprintf (pFile, "\n");

  }


  fprintf (pFile, "%d %d %d  \n",l_olaps.size, l_olaps.exact_size, l_olaps.deleted_size);
  fprintf (pFile, "%d %d %d  \n",l_olaps.exact[0].f_read, l_olaps.exact[0].t_read, l_olaps.exact[0].ori_t);
  
  for (int i=0;i<l_olaps.size;i++)
    fprintf (pFile, "%d %d %d %d %d %d \n",l_olaps.list[i].f_read, l_olaps.list[i].t_read, l_olaps.list[i].ori_t,l_olaps.list[i].offset,l_olaps.list[i].r_arrow,l_olaps.list[i].deleted);
  
  for (int i=0;i<num_of_reads+1;i++)
    fprintf (pFile, "%d \n",l_olaps.f_read_loc[i]);

     
  /*
  for (int i=0;i<list_of_olaps.size;i++){
    fprintf (pFile, " %03d  ",list_of_olaps.list[i].f_read+1);
    fprintf (pFile, "%03d  ",list_of_olaps.list[i].t_read+1);
    if (list_of_olaps.list[i].ori_t)
      fprintf (pFile, "F  ");
    else
      fprintf (pFile, "R  ");
    fprintf (pFile, "%*d\n",4,list_of_olaps.list[i].offset);
    }*/

  fclose (pFile);
#endif


}


void print_unis(unitigs &unis){

}


void find_and_print_unitigs(unitigs unis){
  find_unis(unis);
  print_unis(unis);
};




//
//HW2 codes for finding unitigs end here
//


int main(){

  //HW1 finds overlaps
  read_raw list_of_reads[num_of_reads];
  olaps list_of_olaps;  
  find_and_print_olaps(list_of_reads, list_of_olaps);
  //HW1 ends

  //HW2 finds unitigs
  unitigs unis;
  find_and_print_unitigs(unis);
  //HW2 ends
 
  //HW3 finds a contig
  //HW3 ends


  return 0;
}




/*
 
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
*/

/*
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



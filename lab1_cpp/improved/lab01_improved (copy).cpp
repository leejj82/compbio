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
  int offset;//how much the node is ahead of/behind of the other node.

};

void rc_node(node &a){//reverse complement a read
  a.ori=1-a.ori;
}

class edges_node {
public:
  vector<node> f_node;//nodes from which edges come to this node  
  vector<node> t_node;//nodes to which edges go from this node

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
  int exact_count;//temporary variable for counting used exact overlaps
  
  vector<olap_2> list;//overlap list
  int f_read_loc[num_of_reads+1];//starting point of the f_read in the list_of_olaps
  vector<olap_2> exact;//save the exact overlaps in the list_of_olaps
  edges_node node[num_of_reads];
  
  olaps_2();
};

olaps_2::olaps_2(){
  size=0; deleted_size=0; exact_size=0;exact_size=0;
  fill(f_read_loc, f_read_loc+num_of_reads+1,0);
}

class unitig{
public:
  int size;//# of nodes
  int length;//# length
  vector<node> nodes;//unitig nodes  
  vector<node> f_nodes;//nodes from which the unitig is connected
  vector<node> t_nodes;//nodes to which the unitig is connected
  int f_nodes_size;//# of f_nodes
  int t_nodes_size;//# of t_nodes
  unitig();
};

unitig::unitig(){
  size=0;
  length=read_len;
  f_nodes_size=0;t_nodes_size=0;
}

class unitigs{
public:
  int size;// # of unitigs
  vector<unitig> uni;
  unitigs();
};

unitigs::unitigs(){
  size=0;
}

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
	  temp_node.offset=l_olaps.list[j].offset;
	  l_olaps.node[i].t_node.push_back(temp_node); 

	  l_olaps.node[i].total_edge_ct+=1;
	  l_olaps.node[i].t_edge_ct+=1;

	  if (l_olaps.list[j].ori_t){//second read is front oriented

	    temp_node.num=i;
	    temp_node.ori=1;
	    temp_node.offset=l_olaps.list[j].offset;
	    l_olaps.node[l_olaps.list[j].t_read].f_node.push_back(temp_node);  
	   
	    l_olaps.node[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.node[l_olaps.list[j].t_read].f_edge_ct+=1;	    
	  }
	  else{//second read is reverse complement oriented
	 
	    temp_node.num=i;
	    temp_node.ori=0;
	    temp_node.offset=l_olaps.list[j].offset;
	    l_olaps.node[l_olaps.list[j].t_read].t_node.push_back(temp_node); 
    	
	    l_olaps.node[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.node[l_olaps.list[j].t_read].t_edge_ct+=1;	    
	  }
	}
	else { //location of first read back
	
	  temp_node.num=l_olaps.list[j].t_read;
	  temp_node.ori=l_olaps.list[j].ori_t;
	  temp_node.offset=-l_olaps.list[j].offset;
	  l_olaps.node[i].f_node.push_back(temp_node); 
	
	  l_olaps.node[i].total_edge_ct+=1;
	  l_olaps.node[i].f_edge_ct+=1;

	  if (l_olaps.list[j].ori_t){//second read is front oriented
 
	    temp_node.num=i;
	    temp_node.ori=1;
	    temp_node.offset=-l_olaps.list[j].offset;
	    l_olaps.node[l_olaps.list[j].t_read].t_node.push_back(temp_node); 

	    l_olaps.node[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.node[l_olaps.list[j].t_read].t_edge_ct+=1;	    
	  }
	    
	  
	  else{//second read is reverse complement oriented
	  	   
	    temp_node.num=i;
	    temp_node.ori=0;
	    temp_node.offset=-l_olaps.list[j].offset;
	    l_olaps.node[l_olaps.list[j].t_read].f_node.push_back(temp_node); 
    	
	    l_olaps.node[l_olaps.list[j].t_read].total_edge_ct+=1;
	    l_olaps.node[l_olaps.list[j].t_read].f_edge_ct+=1;	    
	  }
	}
      }
    }
  }

  for (i=0;i<l_olaps.exact_size;i++){
    l_olaps.node[l_olaps.exact[i].t_read].used=1;
  }
}

void find_prior_nodes(node &current_node,olaps_2 &l_olaps,unitig &uni){

  int starting_point=current_node.num; //number of current node
  bool c_forward=current_node.ori;//if the current node is Forward/Reverse Complement (1/0)
  edges_node c_node=l_olaps.node[starting_point];//current node is c_node
  
  node temp_node;
  bool condition;
  
  uni.nodes.push_back(current_node);//insert current node in the unitig
  uni.size++;
  
  if (c_forward){ //current node is in the forward order  
    if (c_node.f_edge_ct>0){ //the current node has incoming edges    
      vector<node> previous_nodes=c_node.f_node;

      condition=(previous_nodes[0].ori  && l_olaps.node[previous_nodes[0].num].t_edge_ct==1)
	||(!previous_nodes[0].ori  && l_olaps.node[previous_nodes[0].num].f_edge_ct==1);//1 if the previous node has one outgoing edge

      if (c_node.f_edge_ct==1 && condition ){ //the current node has exactly one incoming edge and the previous edge has exactly one outgoing edge
	uni.nodes[uni.size-1].offset=previous_nodes[0].offset;
	l_olaps.node[previous_nodes[0].num].used=1;
	find_prior_nodes(previous_nodes[0],l_olaps, uni);
      }
      else{
	uni.f_nodes=previous_nodes;
	uni.f_nodes_size=previous_nodes.size();
      }
    }
  }
  else{ //current node is in the reverse complement order
    if (c_node.t_edge_ct>0){ //the current node has incoming edges    
      vector<node> previous_nodes=c_node.t_node;

      condition=(!previous_nodes[0].ori  && l_olaps.node[previous_nodes[0].num].t_edge_ct==1)
	||(previous_nodes[0].ori  && l_olaps.node[previous_nodes[0].num].f_edge_ct==1);//1 if the previous node has one outgoing edge

      if (c_node.t_edge_ct==1 && condition ){ //the current node has exactly one incoming edge and the previous node has exactly one outgoing edge
	uni.nodes[uni.size-1].offset=previous_nodes[0].offset;
	l_olaps.node[previous_nodes[0].num].used=1;
	node p_node=previous_nodes[0];
	p_node.ori=1-p_node.ori;
	find_prior_nodes(p_node,l_olaps, uni);
      }
      else{
	uni.f_nodes=previous_nodes;
	uni.f_nodes_size=previous_nodes.size();
      }
    }
  }
}

void find_posterior_nodes(node &current_node,olaps_2 &l_olaps,unitig &uni){

  int starting_point=current_node.num; //number of current node
  bool c_forward=current_node.ori;//if the current node is Forward/Reverse Complement (1/0)
  edges_node c_node=l_olaps.node[starting_point];//current node is c_node
  
  node temp_node;
  bool condition;
  
  uni.nodes.push_back(current_node);//insert current node in the unitig
  uni.size++;
  
  if (c_forward){ //current node is in the forward order  
    if (c_node.t_edge_ct>0){ //the current node has outgoing edges    
      vector<node> next_nodes=c_node.t_node;
      
      condition=(next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].f_edge_ct==1)||(!next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].t_edge_ct==1);//1 if the next node has one incoming edge

      if (c_node.t_edge_ct==1 && condition ){ //the current node has exactly one outgoing edge and the next edge has exactly one incoming edge
	l_olaps.node[next_nodes[0].num].used=1;
	find_posterior_nodes(next_nodes[0],l_olaps, uni);
      }
      else{
	uni.t_nodes=next_nodes;
	uni.t_nodes_size=next_nodes.size();
      }
    }
  }
  else{ //current node is in the reverse complement order
    
    if (c_node.f_edge_ct>0){ //the current node has outgoing edges    
      vector<node> next_nodes=c_node.f_node;


      condition=(!next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].f_edge_ct==1)||(next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].t_edge_ct==1);//1 if the next node has one incoming edge
      
      if (c_node.f_edge_ct==1 && condition ){ //the current node has exactly one outgoing edge and the next node has exactly one incoming edge

	l_olaps.node[next_nodes[0].num].used=1;
	node n_node=next_nodes[0];
	n_node.ori=1-n_node.ori;
	find_posterior_nodes(n_node,l_olaps, uni);
      }
      else{
	uni.t_nodes=next_nodes;
	uni.t_nodes_size=next_nodes.size();
      }
    }
  }
} 

void find_a_unitig(int &starting_point, olaps_2 &l_olaps, unitigs &unis){
  
  int i,j;

  node current_node;
  bool condition;
  
  edges_node c_node=l_olaps.node[starting_point];//current node is c_node
  
  if (c_node.total_edge_ct>0){ //there exists at least one edge connected to the node   

    unitig unit;//create a unitig
    unis.uni.push_back(unit); 

    current_node.num=starting_point;
    current_node.ori=1;
    find_prior_nodes(current_node,l_olaps, unis.uni[unis.size]);

    reverse(unis.uni[unis.size].nodes.begin(),unis.uni[unis.size].nodes.end());//reorganize unitig

    if (c_node.t_edge_ct>0){
      vector<node> next_nodes=c_node.t_node;  //consider the next node in the same way   
  
      condition=(next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].f_edge_ct==1)||(!next_nodes[0].ori  && l_olaps.node[next_nodes[0].num].t_edge_ct==1);

      if (c_node.t_edge_ct==1 && condition ){
	
	l_olaps.node[next_nodes[0].num].used=1;
	find_posterior_nodes(next_nodes[0],l_olaps,unis.uni[unis.size]);
      }
      unis.uni[unis.size].t_nodes=next_nodes;
      unis.uni[unis.size].t_nodes_size=next_nodes.size();
    }

    unis.uni[unis.size].nodes[0].offset=0;//set the unitig start from zero

    for (i=0;i<unis.uni[unis.size].size;i++){
      unis.uni[unis.size].length+=unis.uni[unis.size].nodes[i].offset;
    }
    unis.size++;//count the number of unitigs    
    
    /*
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
    */
  }

  else {//single node case
    cout<<"There exists a single node unitig. Not implemented yet.";
  }  
}

void find_unitigs(olaps_2 &l_olaps, unitigs &unis){
  
  int read_exact_match_count=0;
  
  for(int starting_point=0;starting_point<num_of_reads;starting_point++){
    if (!l_olaps.node[starting_point].used){//node not used
      l_olaps.node[starting_point].used=1;
      find_a_unitig(starting_point, l_olaps, unis);
    }
  }
}   

void find_unis(unitigs &unis){

  olaps_2 l_olaps; 
  read_from_olaps(l_olaps);
  record_edge_to_delete(l_olaps);
  set_up_viable_edges(l_olaps);
  // set_up_edges_RC(edges_for_nodes);
  find_unitigs(l_olaps, unis);
 

#if 1
  FILE * pFile;
  pFile = fopen ("lab01.temp","w");

  int s_um=0;
  for(int i=0;i<unis.size;i++){
    s_um+=unis.uni[i].size;
  }
  

  fprintf (pFile, "%d %d \n",unis.size, s_um);  

  for (int i=0;i<unis.size;i++){

    fprintf (pFile, "%d %d %d %d \n", unis.uni[i].size,  unis.uni[i].length,  unis.uni[i].f_nodes_size, unis.uni[i].t_nodes_size );

    fprintf (pFile, "\n");
  }
  
  for (int i=0;i<unis.size;i++){
    for (int j=0;j<unis.uni[i].size;j++){
      fprintf (pFile, "%d %d %d \n", unis.uni[i].nodes[j].num, unis.uni[i].nodes[j].ori, unis.uni[i].nodes[j].offset );
    }
    fprintf (pFile, "\n");
  }

   /*
 
  for (int i=0;i<num_of_reads;i++){
    fprintf (pFile, "%d %d %d %d ",i,l_olaps.node[i].total_edge_ct,l_olaps.node[i].f_edge_ct,l_olaps.node[i].t_edge_ct);
    for (int j=0;j<l_olaps.node[i].f_edge_ct;j++)
      fprintf (pFile, "%d %d %d 1 %d ",  l_olaps.node[i].f_node[j].num,l_olaps.node[i].f_node[j].ori ,i,l_olaps.node[i].f_node[j].offset);
    fprintf (pFile, "\n");
    fprintf (pFile, "%d %d %d %d ",i,l_olaps.node[i].total_edge_ct,l_olaps.node[i].f_edge_ct,l_olaps.node[i].t_edge_ct);

    for (int j=0;j<l_olaps.node[i].t_edge_ct;j++)
      fprintf (pFile, "%d 1 %d %d %d ",i, l_olaps.node[i].t_node[j].num,l_olaps.node[i].t_node[j].ori, l_olaps.node[i].t_node[j].offset);
    fprintf (pFile, "\n");
  }

  fprintf (pFile, "%d %d %d  \n",l_olaps.size, l_olaps.exact_size, l_olaps.deleted_size);
  fprintf (pFile, "%d %d %d  \n",l_olaps.exact[0].f_read, l_olaps.exact[0].t_read, l_olaps.exact[0].ori_t);
  
  for (int i=0;i<l_olaps.size;i++)
    fprintf (pFile, "%d %d %d %d %d %d \n",l_olaps.list[i].f_read, l_olaps.list[i].t_read, l_olaps.list[i].ori_t,l_olaps.list[i].offset,l_olaps.list[i].r_arrow,l_olaps.list[i].deleted);
  
  for (int i=0;i<num_of_reads+1;i++)
    fprintf (pFile, "%d \n",l_olaps.f_read_loc[i]);

     
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

  int i,j;
  FILE * pFile;

#if SAMPLE
  pFile = fopen ("sample.unis","w");
#else
  pFile = fopen ("lab01.unis","w");
#endif
 	
	
  for (i=0;i<unis.size;i++){
 
    fprintf (pFile, "UNI  %02d %*d %*d\n", i+1, 5,unis.uni[i].size,6, unis.uni[i].length);//print the title

    for (j=0;j<unis.uni[i].size;j++){//print the rest
      
      fprintf (pFile, "  %03d  ",unis.uni[i].nodes[j].num+1);
      if (unis.uni[i].nodes[j].ori)
	fprintf (pFile, "F  ");
      else
	fprintf (pFile, "R  ");
      fprintf (pFile, "%*d\n",4,unis.uni[i].nodes[j].offset);
    }
  }
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



  


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
const int l_bd_mp= 700;//950-250
const int u_bd_mp= 760;//1010-250
#else
const int num_of_reads = 300;
const int read_len = 500;
const int l_bd_mp= 1900;//2400-500
const int u_bd_mp= 3100;//3600-500
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

void find_olaps_really(read_raw list_of_reads[num_of_reads], olaps &list_of_olaps){

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

void find_olaps(read_raw list_of_reads[num_of_reads], olaps &list_of_olaps){

  read_from_fasta(list_of_reads); 
  reverse_complement(list_of_reads);
  find_olaps_really(list_of_reads, list_of_olaps);
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

node rc_n(node a){//reverse complement a node
  node temp=a;
  temp.ori=!temp.ori;
  return temp;
}

class edges_node {
public:
  vector<node> f_node;//nodes from which edges come to this node  
  vector<node> t_node;//nodes to which edges go from this node
  
  int total_edge_ct;// #of total edges connected to the node
  int f_edge_ct;// # of edges coming to the node
  int t_edge_ct;//# of edges going out from the node
  bool used;//used=1 if the node is used

  int exact_copy_count;
  edges_node();
};

edges_node::edges_node(){
  total_edge_ct=0;
  f_edge_ct=0;
  t_edge_ct=0;
  used=0;
  exact_copy_count=0;
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

  int f_size;//# of f_nodes/f_unis
  int t_size;//# of t_nodes/t_unis

  unitig();
};

unitig::unitig(){
  size=0;
  length=read_len;
  f_size=0;t_size=0;
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

void edge_setup(int num, bool ori, int offset, edges_node &edge,bool f){

  node temp_node;

  temp_node.num=num;
  temp_node.ori=ori;
  temp_node.offset=offset;
  if (f){
    edge.f_node.push_back(temp_node); 
    edge.f_edge_ct+=1;
  }
  else{
    edge.t_node.push_back(temp_node); 
    edge.t_edge_ct+=1;
  }
  edge.total_edge_ct+=1;
}

void set_up_viable_edges(olaps_2 &l_olaps){

  int i,j;
  
  for (i=0;i<num_of_reads;i++){
    for (j=l_olaps.f_read_loc[i];j<l_olaps.f_read_loc[i+1];j++){
      if (l_olaps.list[j].deleted==0){ //not deleted
	if (l_olaps.list[j].r_arrow){ //location of first read front
	  
	  edge_setup(l_olaps.list[j].t_read,l_olaps.list[j].ori_t, l_olaps.list[j].offset,l_olaps.node[i],0);

	  if (l_olaps.list[j].ori_t)//second read is front oriented
	    edge_setup(i, 1, l_olaps.list[j].offset, l_olaps.node[l_olaps.list[j].t_read],1);	    
	  else //second read is reverse complement oriented
	    edge_setup(i, 0, l_olaps.list[j].offset, l_olaps.node[l_olaps.list[j].t_read],0);	    
	}
	else { //location of first read back

	  edge_setup(l_olaps.list[j].t_read, l_olaps.list[j].ori_t, -l_olaps.list[j].offset, l_olaps.node[i],1);

	  if (l_olaps.list[j].ori_t)//second read is front oriented
	    edge_setup(i,1,-l_olaps.list[j].offset, l_olaps.node[l_olaps.list[j].t_read],0);	    
	  else//second read is reverse complement oriented
	    edge_setup(i, 0, -l_olaps.list[j].offset, l_olaps.node[l_olaps.list[j].t_read],1);	    	  
	}
      }
    }
  }

  for (i=0;i<l_olaps.exact_size;i++){
    l_olaps.node[l_olaps.exact[i].t_read].used=1;
    if(!l_olaps.node[l_olaps.exact[i].f_read].used)
      l_olaps.node[l_olaps.exact[i].f_read].exact_copy_count+=1;
    else
      cout<<" there are multiple exact reads which need to be treated";
  }
}

void insert_node_in_uni (node &current_node, olaps_2 &l_olaps,unitig &uni ){

  l_olaps.node[current_node.num].used=1;//change current_node to used
  uni.nodes.push_back(current_node);//insert current node in the unitig
  uni.size++;

  if (l_olaps.node[current_node.num].exact_copy_count>0){
    node temp_node;
    for (int i=0;i<l_olaps.exact_size;i++){
      temp_node.num=l_olaps.exact[i].t_read;
      if (current_node.ori)
	temp_node.ori=l_olaps.exact[i].ori_t;
      else
	temp_node.ori=1-l_olaps.exact[i].ori_t;
      temp_node.offset=0;
      uni.nodes.push_back(temp_node);
      uni.size++;
    }
  }
}

void insert_boundary_nodes(vector<node> &next_nodes,unitig &uni, bool prior, bool F){
  int i;
  bool RC=!F;
  
  if(prior){//previous nodes
    uni.f_nodes=next_nodes;
    uni.f_size=next_nodes.size();
    if(RC){//reverse complement case
      for (i=0;i<uni.f_size;i++)
	uni.f_nodes[i].ori=!uni.f_nodes[i].ori;
    }
  }
  else{
    uni.t_nodes=next_nodes;
    uni.t_size=(int)(next_nodes.size());
    if(RC){
      for (i=0;i<uni.t_size;i++)
	uni.t_nodes[i].ori=!uni.t_nodes[i].ori;
    }
  }
}


bool one_and_one(edges_node *c_node, olaps_2 &l_olaps,bool prior, bool F){
  bool current_ct;//1 if current node has one edge toward previous or next node
  bool next_ct;//1 if the previous or next node has one edge toward current node
  bool RC=!F;//F=1 current node forward, RC=1 current node RC
  bool forward; //orientation of the previous or next node
  bool t_edge_ct, f_edge_ct;
  
  if(prior && F){
    forward=(*c_node).f_node[0].ori;
    t_edge_ct=l_olaps.node[(*c_node).f_node[0].num].t_edge_ct==1;
    f_edge_ct=l_olaps.node[(*c_node).f_node[0].num].f_edge_ct==1;

    current_ct=(*c_node).f_edge_ct==1;
  }
  else if (prior && RC){
    forward=!(*c_node).t_node[0].ori;
    t_edge_ct=l_olaps.node[(*c_node).t_node[0].num].t_edge_ct==1;
    f_edge_ct=l_olaps.node[(*c_node).t_node[0].num].f_edge_ct==1;
    
    current_ct=(*c_node).t_edge_ct==1;
  }
  else if (!prior && F){
    forward=(*c_node).t_node[0].ori;
    t_edge_ct=l_olaps.node[(*c_node).t_node[0].num].f_edge_ct==1;
    f_edge_ct=l_olaps.node[(*c_node).t_node[0].num].t_edge_ct==1;

    current_ct=(*c_node).t_edge_ct==1;
  }
  else {
    forward=!(*c_node).f_node[0].ori;
    t_edge_ct=l_olaps.node[(*c_node).f_node[0].num].f_edge_ct==1;
    f_edge_ct=l_olaps.node[(*c_node).f_node[0].num].t_edge_ct==1;

    current_ct=(*c_node).f_edge_ct==1;
  }

  next_ct=(forward &&t_edge_ct)||(!forward  && f_edge_ct );
   
  return (current_ct && next_ct);//returns if the previous node has one outgoing edge and the current node has one incoming edge
}

void find_nodes(node &current_node,olaps_2 &l_olaps,unitig &uni,bool prior){

  int i;
  
  int starting_point=current_node.num; //number of current node
  bool c_forward=current_node.ori;//if the current node is Forward/Reverse Complement (1/0)
  edges_node *c_node=&l_olaps.node[starting_point];//current node is c_node
  bool F=1,RC=0; //F: forward node RC: reverse complement node

  if(prior){//previous nodes case 
    if (c_forward){ //current node is in the forward order  
      if ((*c_node).f_edge_ct>0){ //the current node has incoming edges    
	if (one_and_one(c_node, l_olaps,prior,F)){ //the current node has exactly one incoming edge and the previous node has exactly one outgoing edge
	  uni.nodes[uni.size-1].offset=(*c_node).f_node[0].offset;
	  insert_node_in_uni((*c_node).f_node[0], l_olaps, uni);
	  find_nodes((*c_node).f_node[0],l_olaps, uni, prior);
	}
	else
	  insert_boundary_nodes((*c_node).f_node,uni,prior,F);
      }
    }
    else{ //current node is in the reverse complement order
      if ((*c_node).t_edge_ct>0){ //the current node has incoming edges    
	if (one_and_one(c_node, l_olaps,prior,RC)){ //the current node has exactly one incoming edge and the previous edge has exactly one outgoing edge
	  uni.nodes[uni.size-1].offset=(*c_node).t_node[0].offset;
	  node p_node=rc_n((*c_node).t_node[0]);//reverse complement a node
	  insert_node_in_uni(p_node, l_olaps, uni);
	  find_nodes(p_node,l_olaps, uni, prior);
	}
	else
	  insert_boundary_nodes((*c_node).t_node,uni,prior,RC);
      }
    }
  }
  else{//next nodes case
    if (c_forward){ //current node is in the forward order  
      if ((*c_node).t_edge_ct>0){ //the current node has outgoing edges    
	if (one_and_one(c_node, l_olaps,prior,F)){ //the current node has exactly one outgoing edge and the next edge has exactly one incoming edge
	  insert_node_in_uni((*c_node).t_node[0], l_olaps, uni);
	  find_nodes((*c_node).t_node[0],l_olaps, uni, prior);
	}
	else
	  insert_boundary_nodes((*c_node).t_node,uni,prior,F);
      }
    }
    else{ //current node is in the reverse complement order   
      if ((*c_node).f_edge_ct>0){ //the current node has outgoing edges    
      	if (one_and_one(c_node, l_olaps,prior,RC)){ //the current node has exactly one outgoing edge and the next node has exactly one incoming edge
	  node n_node=rc_n((*c_node).f_node[0]);
	  insert_node_in_uni(n_node, l_olaps, uni);
	  find_nodes(n_node,l_olaps, uni, prior);
	}
	else
	  insert_boundary_nodes((*c_node).f_node,uni,prior,RC);
      }
    }
  }
}

void insert_current_node_and_find_prior_and_posterior_nodes (node &current_node,olaps_2 &l_olaps,unitig &uni){

  bool prior=1, posterior=0;
  insert_node_in_uni(current_node, l_olaps, uni);
  find_nodes(current_node,l_olaps, uni, prior);//find the nodes before and including the current node
  reverse(uni.nodes.begin(),uni.nodes.end());//reorder the unitig 
  find_nodes(current_node,l_olaps,uni, posterior);//find the next nodes

  uni.nodes[0].offset=0;//set the unitig starting from zero
  for (int i=0;i<uni.size;i++)//compute unitig lengths
    uni.length+=uni.nodes[i].offset;
}


void find_a_unitig(int &starting_point, olaps_2 &l_olaps, unitigs &unis){
  
  bool condition;
  bool connected_edges_to_the_current_node=l_olaps.node[starting_point].total_edge_ct>0;
 
  if (connected_edges_to_the_current_node){ //there exists at least one edge connected to the node   

    unitig unit;//create a unitig
    unis.uni.push_back(unit); 

    node current_node;
    current_node.num=starting_point;
    current_node.ori=1;

    insert_current_node_and_find_prior_and_posterior_nodes(current_node, l_olaps, unis.uni[unis.size]);
    
    unis.size++;//count the number of unitigs   
  }

  else {//single node case
    cout<<"There exists a single node unitig. Not implemented yet.";
  }  
}

void find_unitigs(olaps_2 &l_olaps, unitigs &unis){
  
  int read_exact_match_count=0;
  
  for(int starting_point=0;starting_point<num_of_reads;starting_point++){
    if (!l_olaps.node[starting_point].used){//node not used
      find_a_unitig(starting_point, l_olaps, unis);
    }
  }
}   

void find_unis(unitigs &unis){

  olaps_2 l_olaps; 
  read_from_olaps(l_olaps);
  record_edge_to_delete(l_olaps);
  set_up_viable_edges(l_olaps);
  find_unitigs(l_olaps, unis);  
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

void find_and_print_unitigs(unitigs &unis){
  find_unis(unis);
  print_unis(unis);
};

//
//HW2 codes for finding unitigs end here
//



//
//HW3 codes for finding contigs start here
//
class con:public unitigs{
public:
  vector<unitig> uni_rc;// reverse complements of unis
  int mate_count;//count the number of mate_pair
  bool used[num_of_reads];//check if the node is used
  
  int length;
  vector<node> unis_list;//unitigs list in contig
  vector<char> contig;

  con();
  void copy_from_unis(unitigs &);
  void set_up_nodes_locations();
};

con::con(){
  mate_count=0;
  fill(used, used+num_of_reads,0);
  length=0;
}

void con::copy_from_unis(unitigs &unis){
  size=unis.size;
  uni=unis.uni;
}

void find_unis_conn (node &Node, unitigs &contig,bool F){
 
  for (int j=0;j<contig.size;j++){
    if (F){
      if(Node.num==contig.uni[j].nodes[contig.uni[j].size-1].num){
	if(Node.ori==contig.uni[j].nodes[contig.uni[j].size-1].ori){
	  Node.num=j;
	  Node.ori=1;
	  break;
	}
      }
      else if(Node.num==contig.uni[j].nodes[0].num){
	if(Node.ori!=contig.uni[j].nodes[0].ori){
	  Node.num=j;
	  Node.ori=0;
	  break;
	}
      }	
    }
    else {
      if (Node.num==contig.uni[j].nodes[contig.uni[j].size-1].num){
	if(Node.ori!=contig.uni[j].nodes[contig.uni[j].size-1].ori){
	  Node.num=j;
	  Node.ori=0;
	  break;
	}
      }
      else if(Node.num==contig.uni[j].nodes[0].num){
	if(Node.ori==contig.uni[j].nodes[0].ori){
	  Node.num=j;
	  Node.ori=1;
	  break;
	}
      }	
    }
  }
}

void set_up_unis_connections(con &contig){//find connections between unitigs and save the info in contig

  int i,k;
  bool F=1,T=0;//F from_unitig T to_unitig

  for (i=0;i<contig.size;i++){
    for (k=0;k<contig.uni[i].f_size;k++)
      find_unis_conn(contig.uni[i].f_nodes[k],contig,F);
    for (k=0;k<contig.uni[i].t_size;k++)
      find_unis_conn(contig.uni[i].t_nodes[k],contig,T);
  }
}

void set_up_uni_rc_and_node_location(con &contig){ //calculate node locations in the unitigs and find and save reverse complement unitigs

  int i,j,temp0,temp1;

  contig.uni_rc=contig.uni;
  
  for( i=0;i<contig.size;i++){

    contig.uni_rc[i].f_size=contig.uni[i].t_size;//take care of connected unitigs for each unitig in reverse complement unitig
    contig.uni_rc[i].f_nodes=contig.uni[i].t_nodes;
    for (j=0;j<contig.uni_rc[i].f_size;j++){
      contig.uni_rc[i].f_nodes[j].ori=!contig.uni_rc[i].f_nodes[j].ori;
    }
    
    contig.uni_rc[i].t_size=contig.uni[i].f_size;
    contig.uni_rc[i].t_nodes=contig.uni[i].f_nodes;
    for (j=0;j<contig.uni_rc[i].t_size;j++){
      contig.uni_rc[i].t_nodes[j].ori=!contig.uni_rc[i].t_nodes[j].ori;
    }

    reverse(contig.uni_rc[i].nodes.begin(),contig.uni_rc[i].nodes.end());//compute reverse complement unitig
    
    for (j=0;j<contig.uni[i].size;j++){
      contig.uni_rc[i].nodes[j].ori=!contig.uni_rc[i].nodes[j].ori;
    }

    temp0=contig.uni_rc[i].nodes[0].offset;   //compute the offset based on read 0
    contig.uni_rc[i].nodes[0].offset=0;   
    for (j=1;j<contig.uni[i].size;j++){
      contig.uni[i].nodes[j].offset+=contig.uni[i].nodes[j-1].offset;
      temp1=contig.uni_rc[i].nodes[j].offset;
      contig.uni_rc[i].nodes[j].offset=contig.uni_rc[i].nodes[j-1].offset+temp0;
      temp0=temp1;
    }
  }
}

void set_up_contig(unitigs &unis, con &contig){
  contig.copy_from_unis(unis);//copy unis to contig
  set_up_unis_connections(contig);//set up connections between unitigs
  set_up_uni_rc_and_node_location(contig);//calculate node locations in the unitigs and find and save reverse complement unitigs
}

void check_self_pairing(con &contig){

  int mate_pair_distance;
  bool zero_one_way, one_zero_way; //check if two reads are numbered in the consecutive order such as (0,1) or (33,32)-possible mate pair

  for (int i=0;i<contig.size;i++){
    for (int j=0;j<contig.uni[i].size-3;j++){ 
      for (int k=j+3;k<contig.uni[i].size;k++){

	zero_one_way=((contig.uni[i].nodes[j].num==contig.uni[i].nodes[k].num-1) && (contig.uni[i].nodes[k].num%2==1));// 0-1 or 10-11 even-odd pair	
	one_zero_way=((contig.uni[i].nodes[j].num==contig.uni[i].nodes[k].num+1) && (contig.uni[i].nodes[k].num%2==0));//1-0 or 15-14 odd-even pair

	if( zero_one_way || one_zero_way ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if( contig.uni[i].nodes[j].ori==1 && contig.uni[i].nodes[k].ori==0 ){// two reads are facing each other 5'-3' 3'-5' way    

	    mate_pair_distance=contig.uni[i].nodes[k].offset-contig.uni[i].nodes[j].offset;

	    if ((l_bd_mp <=mate_pair_distance ) && (mate_pair_distance <=u_bd_mp )){//two reads are distanced between 2400-3600
	      contig.mate_count+=1;
	      contig.used[contig.uni[i].nodes[j].num]=1;
	      contig.used[contig.uni[i].nodes[k].num]=1;
	    }
	  }
	}
      }
    }
  }
}

void really_find_contig(con &contig){
    
  int i,start_uni=0;

  for (i=0;i<contig.size;i++){//if an end of a unitig does not have connecting edges, then the unitig can be a boundary of a contig
    if(contig.uni[i].f_size==0 || contig.uni[i].t_size==0){
      if(contig.uni[i].f_size+contig.uni[i].t_size==0)
	cout<<"There is a unitig without any edges to outside"<<endl;
      else{
	start_uni=i;
	break;
      }
    }
  }

  node first_uni;
  first_uni.num=start_uni;
  first_uni.ori=1;
  first_uni.offset=0;
  contig.unis_list.push_back(first_uni);//insert the first unitig in contig unitig list
  
  /*  if (iterate_for_finding_a_contig(contig))
    cout<<"found a contig"<<endl;
  else
  cout<<"could not find a contig"<<endl;*/

}

void find_contig(read_raw list_of_reads[num_of_reads], unitigs &unis,con &contig){

  set_up_contig(unis, contig);//copy relevant info from unis, create reverse complement unitigs, find connections between unitigs 
  check_self_pairing(contig);
  really_find_contig(contig);
}

void print_contig(con &contig){ 

}

void find_and_print_a_contig(read_raw list_of_reads[num_of_reads], unitigs &unis,con &contig){
  find_contig(list_of_reads, unis, contig);
  print_contig(contig);

#if 1
  cout<<contig.mate_count;
  FILE * pFile;
  pFile = fopen ("lab01.temp","w");

  for (int i=0;i<unis.size;i++){
    for (int j=0;j<contig.uni[i].size;j++){
      fprintf (pFile, "%d %d %d   ", contig.uni[i].nodes[j].num, contig.uni[i].nodes[j].ori, contig.uni[i].nodes[j].offset );
      fprintf (pFile, "%d %d %d   ", contig.uni_rc[i].nodes[j].num, contig.uni_rc[i].nodes[j].ori, contig.uni_rc[i].nodes[j].offset );

      fprintf (pFile, "\n");
    }
    fprintf (pFile, "\n\n\n\n\n");
  }
 
  for (int i=0;i<unis.size;i++){
    /*    fprintf (pFile, "%d %d \n", unis.uni[i].f_size, unis.uni[i].t_size );
    for (int j=0;j<unis.uni[i].f_size;j++){
      fprintf (pFile, "%d %d \n", unis.uni[i].f_nodes[j].num, unis.uni[i].f_nodes[j].ori );
    }
    for (int j=0;j<unis.uni[i].t_size;j++){
      fprintf (pFile, "%d %d \n", unis.uni[i].t_nodes[j].num, unis.uni[i].t_nodes[j].ori );
    }
    fprintf (pFile, "\n");*/

    fprintf (pFile, "%d %d \n", contig.uni[i].f_size, contig.uni[i].t_size );
    for (int j=0;j<contig.uni[i].f_size;j++){
      fprintf (pFile, "%d %d \n", contig.uni[i].f_nodes[j].num, contig.uni[i].f_nodes[j].ori );
    }
    for (int j=0;j<contig.uni[i].t_size;j++){
      fprintf (pFile, "%d %d \n", contig.uni[i].t_nodes[j].num, contig.uni[i].t_nodes[j].ori );
    }
    fprintf (pFile, "\n");

    
    fprintf (pFile, "%d %d \n", contig.uni_rc[i].f_size, contig.uni_rc[i].t_size );
    for (int j=0;j<contig.uni_rc[i].f_size;j++){
      fprintf (pFile, "%d %d \n", contig.uni_rc[i].f_nodes[j].num, contig.uni_rc[i].f_nodes[j].ori );
    }
    for (int j=0;j<contig.uni_rc[i].t_size;j++){
      fprintf (pFile, "%d %d \n", contig.uni_rc[i].t_nodes[j].num, contig.uni_rc[i].t_nodes[j].ori );
    }
    fprintf (pFile, "\n");
  }

  fprintf (pFile, "\n\n\n\n\n");

  fclose(pFile);
#endif
}



//
//HW3 codes for finding contigs end here
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
  con contig;
  find_and_print_a_contig(list_of_reads,unis,contig);
  //HW3 ends


  return 0;
}

/* 

int mate_pair_check_1(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC, vector<vector<int> > &unitigs_info, int unitigs_con_count[][3], int uni1, int uni1_FR, int uni2, int uni2_FR, int distance){

  int i,j,mate_pair_distance;
  int mate_pair_count=0;
  
  for (i=0;i<unitigs_info[uni1][0];i++){
    for (j=0;j<unitigs_info[uni2][0];j++){

      if (uni1_FR==1 && uni2_FR==1){
	if( ((unis[uni1][i][0]==unis[uni2][j][0]+1) && (unis[uni2][j][0]%2==0)) || ((unis[uni1][i][0]==unis[uni2][j][0]-1) && (unis[uni2][j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	  if(unis[uni1][i][1]==1 && unis[uni2][j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  
	    mate_pair_distance=distance+unitigs_info[uni1][1]-unis[uni1][i][2]+unis[uni2][j][2];

	    if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) ){//two reads are distanced between 2400-3600
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
	      mate_pair_count+=1;
	    }	  
	  }	  
	}
      }
    }
  }
  return mate_pair_count;
}

// Driver function to sort the 2D vector
// on basis of a particular column
bool sortcol( const vector<int>& v1,
	      const vector<int>& v2 ) {
  return v1[3] >  v2[3];
}

int check_for_validity(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info,int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con, vector<vector<int> > &contig_unis_list, int &mate_count, vector<vector<int> > &real_contig_unis_list, vector<vector<int> > &contig_reads){

  int i,j,k;
  int total_num_of_reads_in_contig=0;
  int FB;
  int contig_unis_len=contig_unis_list.size();
  int sum, mate_pair_distance;
  int reads_check[num_of_reads+1]={0};
  int contig_len;
  
  for (i=0;i<contig_unis_len;i++){
    total_num_of_reads_in_contig+=unitigs_info[contig_unis_list[i][0]][0];
  }

  int contig_reads_all[total_num_of_reads_in_contig][4]; //read#, F/RC, starting position, mate_pair yes/no 1/0

  
  int contig_size=0;
  for (i=0;i<contig_unis_len;i++){

    FB=contig_unis_list[i][1];
    if (i==0)
      sum=0;
    else{
      sum=sum+unitigs_info[contig_unis_list[i-1][0]][1]+contig_unis_list[i][2];
    }


    if (FB==1){//front unitig
      for (j=0;j<unitigs_info[contig_unis_list[i][0]][0];j++){
	contig_reads_all[contig_size][0]=unis[contig_unis_list[i][0]][j][0];
	contig_reads_all[contig_size][1]=unis[contig_unis_list[i][0]][j][1];
	contig_reads_all[contig_size][2]=unis[contig_unis_list[i][0]][j][2]+sum;
	contig_reads_all[contig_size][3]=0;
	contig_size++;
      }
    }
    else{ //reverse complement unitig
      for (j=0;j<unitigs_info[contig_unis_list[i][0]][0];j++){
	contig_reads_all[contig_size][0]=unis_RC[contig_unis_list[i][0]][j][0];
	contig_reads_all[contig_size][1]=unis_RC[contig_unis_list[i][0]][j][1];
	contig_reads_all[contig_size][2]=unis_RC[contig_unis_list[i][0]][j][2]+sum;
	contig_reads_all[contig_size][3]=0;	
	contig_size++;
      }
    }
  }  
  contig_len=contig_reads_all[contig_size-1][2];//in fact, contig_length-read_length
  
  
  for (i=0;i<total_num_of_reads_in_contig-1;i++){
    for (j=i+1;j<total_num_of_reads_in_contig;j++){

      if( ((contig_reads_all[i][0]==contig_reads_all[j][0]+1) && (contig_reads_all[j][0]%2==0)) || ((contig_reads_all[i][0]==contig_reads_all[j][0]-1) && (contig_reads_all[j][0]%2==1)) ){//two reads are in the consecutive order such as (0,1) or (33,32)-possible mate pair
	if(contig_reads_all[i][1]==1 && contig_reads_all[j][1]==0 ){// two reads are facing each other 5'-3' 3'-5' way
	  mate_pair_distance=contig_reads_all[j][2]-contig_reads_all[i][2];      
	  if ( (l_bd_mp <=mate_pair_distance) && (u_bd_mp >=mate_pair_distance) &&( contig_reads_all[i][3]==0 && contig_reads_all[j][3]==0) ){//two reads are distanced between 2400-3600 and not used
	    contig_reads_all[i][3]=1;contig_reads_all[j][3]=1;
	    if ( reads_check[contig_reads_all[i][0]]==0 && reads_check[contig_reads_all[j][0]]==0){
	      reads_check[contig_reads_all[i][0]]=1;reads_check[contig_reads_all[j][0]]=1;reads_check[num_of_reads]+=2;
	    }
	  }
	}
      }  
    }
  }

  
  int temp1,temp2;

  for (i=0;i<total_num_of_reads_in_contig;i++){
    if( contig_reads_all[i][3]==0){
      for (j=1;j<i+1;j++){
	if ((contig_reads_all[i-j][3]==1) && ((contig_reads_all[i][2]-contig_reads_all[i-j][2])<read_len) ){
	  temp1=i-j;
	  break;
	}
	else if( (contig_reads_all[i][2]-contig_reads_all[i-j][2])>=read_len){
	  if ( contig_len-contig_reads_all[i][2]>u_bd_mp){
	    return -1;
	  }
	  else{
	    return 0;
	  }
	}
      }
      for (k=i+1;k<total_num_of_reads_in_contig;k++){
	if (contig_reads_all[k][3]==1 && (contig_reads_all[k][2]-contig_reads_all[i][2])<read_len){
	  temp2=k;
	  break;
	}
	else if( (contig_reads_all[k][2]-contig_reads_all[i][2])>=read_len){
	  if ( contig_len-contig_reads_all[i][2]>u_bd_mp){
	    return -1;
	  }
	  else{
	    return 0;
	  }
	}
      }
      if( contig_reads_all[temp2][2]-contig_reads_all[temp1][2]>=read_len){
	if ( contig_len-contig_reads_all[i][2]>u_bd_mp){
	  return -1;
	}
	else{
	  return 0;
	}

      }
    }
  }
 

  
  if (reads_check[num_of_reads]==num_of_reads){
    real_contig_unis_list=contig_unis_list;

    vector<int> temp(3);
    for (i=0;i<total_num_of_reads_in_contig;i++){
      if(contig_reads_all[i][3]!=0){
	for (j=0;j<3;j++){
	  temp[j]=contig_reads_all[i][j];  
	}
	contig_reads.push_back(temp);
      }
    }
    return 1;
  }
  else {
    return 0;
  }
}

int iterate_for_finding_a_contig(vector<vector<vector<int> > > &unis, vector<vector<vector<int> > > &unis_RC,int &num_of_unitigs,vector<vector<int> > &unitigs_info,int unitigs_con_count[][3], vector<vector<vector<vector<int> > > > &unitigs_con, vector<vector<int> > &contig_unis_list, int &mate_count, vector<vector<int> > &real_contig_unis_list, vector<vector<int> > &contig_reads){
  
  int contig_size=(int)(contig_unis_list.size());
  int i,j,k;
  int uni1,uni1_FR,uni2,uni2_FR;
  int last_uni, last_uni_FR;
  
  int distance, distance0;
  int mate_pair_count;
  int result1, result2;
  
  
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
	distance0=-read_len;
      else
	distance0=distance0+unitigs_info[contig_unis_list[j+1][0]][1]+contig_unis_list[j+1][2];
    
      distance=distance0+unitigs_con[last_uni][last_uni_FR][k][4];

      if (distance<=u_bd_mp ){
	mate_pair_count+= mate_pair_check_1(unis, unis_RC, unitigs_info, unitigs_con_count, uni1, uni1_FR, uni2, uni2_FR, distance);
      }
    }

    temp[0]=uni2;
    temp[1]=uni2_FR;
    temp[2]=-read_len+unitigs_con[last_uni][last_uni_FR][k][4];//distance between contig and the next unitig(negative number of course)
    temp[3]=mate_pair_count;
    next_uni_table.push_back(temp);

  }    

  sort(next_uni_table.begin(), next_uni_table.end(),sortcol);
  
  for(i=0;i<next_uni_table.size();i++){
    vector<vector<int> > contig_unis_list_2=contig_unis_list;
    contig_unis_list_2.push_back(next_uni_table[i]);
    int mate_count_2=mate_count+next_uni_table[i][3];
    result1=check_for_validity(unis, unis_RC, num_of_unitigs,unitigs_info,unitigs_con_count, unitigs_con, contig_unis_list_2, mate_count_2, real_contig_unis_list, contig_reads);
    
    if(result1==0){
      result2=iterate_for_finding_a_contig(unis, unis_RC, num_of_unitigs,unitigs_info,unitigs_con_count, unitigs_con, contig_unis_list_2, mate_count_2, real_contig_unis_list, contig_reads);
      if (result2==1){
	return result2;
      }
    }
    else if (result1==1){
      return result1;
    }
  }
  
return 0;
}






void record_and_print_a_contig(char *contig, int extra[][2], vector<vector<int> > &contig_reads, char list_of_reads[][read_len+1], char list_of_reads_RC[][read_len+1]){
  int i;
  int contig_len=contig_reads[num_of_reads-1][2]+read_len;
  
  for (i=0;i<num_of_reads-1;i++){
    if (contig_reads[i][1]==1){
      strncpy(contig+contig_reads[i][2], list_of_reads[contig_reads[i][0]],contig_reads[i+1][2]-contig_reads[i][2]);
    }
    else{ //reverse complement order
      strncpy(contig+contig_reads[i][2], list_of_reads_RC[contig_reads[i][0]],contig_reads[i+1][2]-contig_reads[i][2]);
    }
  }

  if (contig_reads[i][1]==1){
    strncpy(contig+contig_reads[num_of_reads-1][2], list_of_reads[contig_reads[num_of_reads-1][0]],read_len);
  }
  else{ //reverse complement order
    strncpy(contig+contig_reads[num_of_reads-1][2], list_of_reads_RC[contig_reads[num_of_reads-1][0]],read_len);
  }
  
  contig[contig_len]='\0';


  for (i=0;i<num_of_reads;i++){
    if (contig_reads[i][1]==1){
      extra[contig_reads[i][0]][0]=contig_reads[i][2];
      extra[contig_reads[i][0]][1]=contig_reads[i][2]+read_len-1;
    }
    else{//RC case
      extra[contig_reads[i][0]][0]=contig_reads[i][2]+read_len-1;
      extra[contig_reads[i][0]][1]=contig_reads[i][2];
    }
  }

  FILE * pFile;

  pFile = fopen ("lab01.contig","w");
  fprintf (pFile,">Contig\n");
  fprintf (pFile, "%s  ",contig);
  fclose (pFile);

  pFile = fopen ("lab01.extra","w");

  for(i=0;i<num_of_reads;i++){
    fprintf (pFile, "%03d  ",i+1);
    fprintf (pFile, "%*d  ",5,extra[i][0]);
    fprintf (pFile, "%*d\n",5,extra[i][1]);  
}
  fclose (pFile);
 
}

int main(){
  vector<vector<int> > real_contig_unis_list; //unitig #, F/RC, starting pt
  vector<vector<int> > contig_reads;//read#, F/RC, starting pt
  find_a_contig(unitigs,unitigs_info,real_contig_unis_list, contig_reads);
  
  char contig[contig_reads[num_of_reads-1][2]+read_len+1];
  int extra[num_of_reads][2];
  record_and_print_a_contig(contig,extra, contig_reads,list_of_reads,list_of_reads_RC);
	


 */




  /*
#if 1
  FILE * pFile;
  pFile = fopen ("lab01.temp","w");

  

 
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
    }

  fclose (pFile);
  #endif
*/

  


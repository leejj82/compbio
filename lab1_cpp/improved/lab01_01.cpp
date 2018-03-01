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

class READ{
public:
  char read[read_len+1]; //read
  char read_rc[read_len+1]; //reverse complement of the read
  READ();
  void reverse_complement();
};

READ::READ () {
  read[0]='\0';
}

void READ::reverse_complement(){
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
  int from_read,to_read;
  int ori; // orientation=1/0 if to_read is forward/reverse complement
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

void read_from_fasta(READ list_of_reads[num_of_reads]){

  char str[read_len];
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
	    strcat(list_of_reads[i].read,str);
	  }
      }
  }
  infile.close();
}

void reverse_complement(READ list_of_reads[num_of_reads]){
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
      i++;matched_len++;
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
	return 0;
      }
    }
  }
  return 1;
}

bool find_overlaps_of_two_reads(READ &first, READ &second,olap &Olap, read_ends &Read_ends){

  bool not_found=1;

  if (not_found){ 
    not_found=KMP_search(first.read, Read_ends.right, second.read, Olap);
    Olap.ori=1;
    if (not_found){
      not_found=KMP_search(first.read, Read_ends.right, second.read_rc, Olap);
      Olap.ori=0;
      if (not_found){
	not_found=KMP_search(first.read, Read_ends.left, second.read, Olap);
	Olap.ori=1;
	if (not_found){
	  not_found=KMP_search(first.read, Read_ends.left, second.read_rc, Olap);
	  Olap.ori=0;
	}
      }
    }
  }
  return 1-not_found;
}

void find_olaps(READ list_of_reads[num_of_reads], olaps &list_of_olaps){

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
      Olap.from_read=first_read;
      Olap.to_read=second_read;
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
    fprintf (pFile, " %03d  ",list_of_olaps.list[i].from_read+1);
    fprintf (pFile, "%03d  ",list_of_olaps.list[i].to_read+1);
    if (list_of_olaps.list[i].ori==1)
      fprintf (pFile, "F  ");
    else
      fprintf (pFile, "R  ");
    fprintf (pFile, "%*d\n",4,list_of_olaps.list[i].offset);
  }
  
  fclose (pFile);
}

void find_and_print_olaps(READ list_of_reads[num_of_reads], olaps &list_of_olaps){

  read_from_fasta(list_of_reads); 
  reverse_complement(list_of_reads);
  find_olaps(list_of_reads,list_of_olaps);
  print_olaps(list_of_olaps);
}

//
//HW1 codes for finding overlaps end here
//

int main(){


  //HW1 finds overlaps
  READ list_of_reads[num_of_reads];
  olaps list_of_olaps;  
  find_and_print_olaps(list_of_reads,list_of_olaps);
  //HW1 ends
  return 0;
}

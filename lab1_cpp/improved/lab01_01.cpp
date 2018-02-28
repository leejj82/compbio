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

void KMP_table(char *pat,int *failure){

  int cnd;
  failure[0] = -1;

  for (int pos=1;pos < 40;pos++){
    cnd=failure[pos-1];

    while((pat[pos]!=pat[cnd+1]) && (cnd>=0)){
      cnd=failure[cnd];
    }

    if (pat[pos]==pat[cnd+1]){
      failure[pos]=cnd+1;
    }
    else {
      failure[pos]=-1;
    }
  }
}

bool check_rest(int LEFT, int match_start, int match_end, char *list_original, char *list_compare){
  int i;
  if (LEFT){//left
    for(i=match_end+1;i<read_len;i++){
      if (list_original[i-match_start]!=list_compare[i])
	return 0;
    }
  }
  else {//right
    for (i=0;i<match_start;i++){
      if (list_original[read_len-match_end+i-1]!=list_compare[i])
	return 0;
    }
  }
  return 1;
}

bool KMP_search(char *list_original, char *pat, bool LEFT, char *list_compare, int *failure, int *Overlap){

  int i=0,j,matched_len=0, match_start,match_end;

  while (i < read_len){
    if (pat[matched_len] == list_compare[i]){
      i++;matched_len++;
    }
    else if (matched_len==0)
      i++;
    else
      matched_len=failure[matched_len-1]+1;
    if (matched_len==40){
      match_start=i-40;
      match_end=i-1;
      if (check_rest(LEFT, match_start,match_end,list_original,list_compare)){
	*Overlap=-LEFT*match_start+(1-LEFT)*(read_len-match_end-1);
	return 0;
      }
    }
  }
  return 1;
}

bool find_overlaps_of_two_reads(char list_original[read_len + 1], char list_compare[read_len + 1], char list_compare_rc[read_len+1],int *Forward_Backward,int *Overlap, int * failure_R, int * failure_L, char * right_end, char * left_end){

  bool LEFT=1;
  bool not_found=1;

  if (not_found){ 
    not_found=KMP_search(list_original, right_end, !LEFT, list_compare, failure_R, Overlap);
    *Forward_Backward=1;
    if (not_found){
      not_found=KMP_search(list_original, right_end, !LEFT, list_compare_rc, failure_R, Overlap);
       *Forward_Backward=0;
      if (not_found){
	not_found=KMP_search(list_original, left_end, LEFT, list_compare, failure_L, Overlap);
	 *Forward_Backward=1;
	if (not_found){
	  not_found=KMP_search(list_original, left_end, LEFT, list_compare_rc, failure_L, Overlap);
	   *Forward_Backward=0;
	}
      }
    }
  }
  return 1-not_found;
}

void find_olaps(READ list_of_reads[num_of_reads], olaps &list_of_olaps){

  int i, j;
  int Forward_Backward[1], Overlap[1];
  
  char right_end[41];
  char left_end[41];

  int failure_R[40];
  int failure_L[40];

  olap temp;
  
  list_of_olaps.size=0;
  
  for (i=0;i<num_of_reads-1;i++){

    strncpy(right_end, list_of_reads[i].read+read_len-40,40);
    strncpy(left_end, list_of_reads[i].read,40);
    right_end[40]='\0';
    left_end[40]='\0';
    KMP_table(right_end, failure_R);
    KMP_table(left_end, failure_L);    	

    for (j=i+1;j<num_of_reads;j++){
      if (find_overlaps_of_two_reads(list_of_reads[i].read,list_of_reads[j].read,list_of_reads[j].read_rc,Forward_Backward,Overlap,failure_R, failure_L,right_end,left_end)){
	temp.from_read=i;
	temp.to_read=j;
	temp.ori=*Forward_Backward;
	temp.offset=*Overlap;
	list_of_olaps.size++;
	list_of_olaps.list.push_back(temp);
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

int main(){

  READ list_of_reads[num_of_reads];
  olaps list_of_olaps;  
  find_and_print_olaps(list_of_reads,list_of_olaps);

  return 0;
}

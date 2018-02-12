#include<bits/stdc++.h>
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

#if SAMPLE
class Edge
{
public:
  int a;
  int b;
};

struct mycmp {
  bool operator() (const Edge& e1, const Edge& e2) {
    if(e1.a < e2.a)
      return true;
    else if(e2.a < e1.a)
      return false;

    return false;
  }
} myobject;

void test_sort() {
  Edge edges[2];
  edges[0].a = 10;
  edges[1].a = 5;

  sort(edges, edges+2, myobject);

  cout << "1: " << edges[0].a << endl;
  cout << "2: " << edges[1].a << endl;
}
#endif

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

int check_rest(int RIGHT_LEFT, int match_start, int match_end, char *list_original, char *list_compare){
  int i;
  if (RIGHT_LEFT==0){
    for (i=0;i<match_start;i++){
      if (list_original[read_len-match_end+i-1]!=list_compare[i])
	return 0;
    }
  }
  else {
    for(i=match_end+1;i<read_len;i++){
      if (list_original[i-match_start]!=list_compare[i])
	return 0;
    }
  }

  return 1;

}

int KMP_search(char *list_original, char *pat,int RIGHT_LEFT, char *list_compare, int *failure, int *Overlap){

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
      if (check_rest(RIGHT_LEFT, match_start,match_end,list_original,list_compare)==1){
	*Overlap=-RIGHT_LEFT*match_start+(1-RIGHT_LEFT)*(read_len-match_end-1);
	return 1;
      }
    }
  }
  return 0;
}

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

#if DEBUG // || 1
  int num_reads = 0;
  for(int i = 0; i < num_of_reads; i++) {
    char* str = list_of_reads[i];
    int len = strlen(str);
    if(len <= 0) break;
    if(len != 500) {
      cerr << "Error: read " << i + 1 << " is not correct." << endl;
    }
    num_reads++;
  }

  cerr << "Number of reads: " << num_reads << endl;
#endif

  /*
#if 0
  int a = 0;
  int b = 0;
#endif
  */
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

int find_overlaps_of_two_reads(char list_original[read_len + 1], char list_compare[read_len + 1], char list_compare_rc[read_len+1],int *Forward_Backward,int *Overlap){

  char right_end[41];
  char left_end[41];
  int  failure[40];
  int RIGHT=0;
  int LEFT=1;
  int found=0;

  strncpy(right_end, list_original+read_len-40,40);
  strncpy(left_end, list_original,40);
  right_end[40]='\0';
  left_end[40]='\0';
    	
  if (found==0){ 
    KMP_table(right_end, failure);
    found=KMP_search(list_original, right_end, RIGHT, list_compare, failure, Overlap);
    *Forward_Backward=1;
    if (found==0){
      found=KMP_search(list_original, right_end, RIGHT, list_compare_rc, failure, Overlap);
       *Forward_Backward=0;
      if (found==0){
	KMP_table(left_end, failure);
	found=KMP_search(list_original, left_end, LEFT, list_compare, failure, Overlap);
	 *Forward_Backward=1;
	if (found==0){
	  found=KMP_search(list_original, left_end, LEFT, list_compare_rc, failure, Overlap);
	   *Forward_Backward=0;
	}
      }
    }
  }
  return found;
}

int find_overlaps(char list_of_reads[][read_len+1],char list_of_reads_RC[][read_len+1],int list_of_overlaps[][4]){
  int i, j ,k=0;
  int Forward_Backward[1], Overlap[1];
  for (i=0;i<num_of_reads-1;i++){
    for (j=i+1;j<num_of_reads;j++){
      if (find_overlaps_of_two_reads(list_of_reads[i],list_of_reads[j],list_of_reads_RC[j],Forward_Backward,Overlap)==1){
	list_of_overlaps[k][0]=i;
	list_of_overlaps[k][1]=j;
	list_of_overlaps[k][2]=*Forward_Backward;
	list_of_overlaps[k][3]=*Overlap;
	k++;
      }
    }
  }
  return k;
}

void print_overlaps(int number_of_overlaps, int list_of_overlaps[][4]){
  int i,ind;

  FILE * pFile;

  #if SAMPLE
  pFile = fopen ("sample.olaps","w");
  #else
  pFile = fopen ("lab1.olaps","w");
  #endif

  for (i=0;i<number_of_overlaps;i++){
    fprintf (pFile, " %03d  ",list_of_overlaps[i][0]+1);
    fprintf (pFile, "%03d  ",list_of_overlaps[i][1]+1);
    ind=list_of_overlaps[i][2];  
    if (ind==1)
      fprintf (pFile, "F  ");
    else
      fprintf (pFile, "R  ");
    fprintf (pFile, "%*d\n",4,list_of_overlaps[i][3]);
  }
  
  fclose (pFile);
}


int main()
{
  int i, j;

  //define reads and get an input from a file
  //int num_of_reads=300;

  #if SAMPLE
  test_sort();
  #endif

  char list_of_reads[num_of_reads][read_len + 1];
  for (i=0;i<num_of_reads;i++){
      // list_of_reads[i]=(char*)malloc(501);
      list_of_reads[i][0] = '\0';
    }

  read_from_fasta(list_of_reads);

  int read_length=strlen(list_of_reads[0]);

  //define reverse complement reads
  char  list_of_reads_RC[num_of_reads][read_len + 1];
  for (i=0;i<num_of_reads;i++){
      list_of_reads_RC[i][0] = '\0';
    }
  
  reversecomplement(list_of_reads, list_of_reads_RC);

  //find overlaps and save it to a list
  int list_of_overlaps[num_of_reads*(num_of_reads-1)/2][4];
  
  int number_of_overlaps=find_overlaps(list_of_reads,list_of_reads_RC, list_of_overlaps);

  print_overlaps(number_of_overlaps,list_of_overlaps);

  cout<<number_of_overlaps<<"\n";

  return 0;
}



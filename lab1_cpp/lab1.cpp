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

void KMP_table(char *pat, int M, int *failure)
{
  int pos = 1;
  int cnd = 0;
  failure[0] = -1;

  while (pos < M)
  {
    if (pat[pos] == pat[cnd])
    {
      failure[pos]=failure[cnd];
      pos++;
      cnd++;
    }
    else
    {
      failure[pos]=cnd;
      cnd=failure[cnd];
      while (cnd >=0 and pat[pos]!=pat[cnd])
      {
	cnd=pat[cnd];
      }
      pos++;
      cnd++;
    }
  }
  failure[pos]=cnd;
}

void KMP_search(char *pat, char *txt)
{
  int M = strlen(pat);
  int N = strlen(txt);
  int failure[M];

  KMP_table(pat, M, failure);

  int m = 0;
  int i = 0;
  while (m+i < N)
  {
    if (pat[i] == txt[m+i])
    {
      if (i+1 == M)
      {
	cout<<"Found pattern at index "<<m<<endl;
	for (int k=m;k<m+M;k++)
	  {
	    cout<<txt[k];
	  }
	cout<<endl;
	for (int k=0;k<M;k++)
	  {
	    cout<<failure[k];
	  }
	cout<<endl;
      	m=m+i-failure[i];
	i=failure[i];
      }
      else
      {
	i=i+1;
      }
    }
    else
    {
      if (failure[i]>-1)
	{
	  m=m+i-failure[i];
	  i=failure[i];	  
	}
      else
	{
	  m=m+i+1;
	  i=0;
	}
    }
  }
}

int pattern_example()
{
  char *txt = "ABABDABACDABABCABABDSEABABCABAB";
  char *pat = "ABABCABAB";
  KMP_search(pat, txt);
  return 0;
}


void read_from_fasta(char list_of_reads[][read_len + 1])
{
  char str[500];
  //char *str= temp;
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

#if DEBUG || 1
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

#if 0
  int a = 0;
  int b = 0;
#endif
}

void reversecomplement(int num_of_reads,int read_length, char list[][read_len + 1], char list_rc[][read_len + 1]){
  for (int j=0;j<num_of_reads;j++)
    {
      for (int i=0;i<read_length;i++)
	{
	  if (list[j][read_length-i-1]=='a')
	    {
	      list_rc[j][i]='t';
	    }
	  else if (list[j][read_length-i-1]=='t')
	    {
	      list_rc[j][i]='a'; 
	    }
	  else if (list[j][read_length-i-1]=='c')
	    {
	      list_rc[j][i]='g';
	    }
	  else if (list[j][read_length-i-1]=='g')
	    {
	      list_rc[j][i]='c';
	    }
	}
      list_rc[j][read_length]='\0';
    }
 }

/*
void find_overlaps(int num_of_reads,int read_length, char **list, char **list_rc)
{
  for (int i=0;i<num_of_reads;i++)
    {
      for (int j=i+1;j<num_of_reads;j++)
	{
	  find_overlaps_of_two_reads(read_length, list[i],list[j],list_rc[j]);
	}
    }
}
*/

void find_overlaps_of_two_reads(int read_length, char *list_original, char *list_compare, char *list_compare_rc)
{
  char right_end[21];
  char left_end[21];
  for (int i=0;i<read_length;i++)
    {
      strncpy(right_end, list_original+read_length-20,20);
      strncpy(left_end, list_original,20);
      right_end[20]='\0';
      left_end[20]='\0';
      cout<<list_original<<"\n"<<right_end<<strlen(right_end)<<"\n"<<left_end<<strlen(left_end)<<"\n";
    }
}
	  
int main()
{
  //define reads and get an input from a file
  //int num_of_reads=300;

  test_sort();
  exit(1);
  
  char list_of_reads[num_of_reads][read_len + 1];
  for (int i=0;i<num_of_reads;i++)
    {
      // list_of_reads[i]=(char*)malloc(501);
      list_of_reads[i][0] = '\0';
    }

  read_from_fasta(list_of_reads);

  int read_length=strlen(list_of_reads[0]);

  //define reverse complement reads
  char  list_of_reads_RC[num_of_reads][read_len + 1];
  for (int i=0;i<num_of_reads;i++)
    {
      list_of_reads_RC[i][0] = '\0';
    }
  
  reversecomplement(num_of_reads, read_length, list_of_reads, list_of_reads_RC);

  /*  //find overlaps
  find_overlaps(num_of_reads,read_length,list_of_reads,list_of_reads_RC);
  */
  find_overlaps_of_two_reads(read_length, list_of_reads[0],list_of_reads[1],list_of_reads_RC[1]);
    
  return 0;
}



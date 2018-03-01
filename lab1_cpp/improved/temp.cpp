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
  edge_b r();
  edge_b rc();

};

edge_b edge_b::r(){
  edge_b temp;
  temp.from_read=to_read;
  temp.to_read=from_read;
  return temp;
}
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

int main(){ 

  return 0;
}

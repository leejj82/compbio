//#include<bits/stdc++.h>
#include<stdlib.h>
#include<cstring>
#include<algorithm>

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class read{
public:
  int num[2];//=read number 0-299
};


int main(){ 

  read a,b,c;
  a.num[0]=3;a.num[1]=1;
  b=a;
  b.num[0]=2;

  cout<<a.num[0];
  return 0;
}

#include<vector>
#include<iostream>

using namespace std;

int main(){
 
  vector<int> v;

  v.reserve(10);

  v[0]=4;
  v.push_back(3);

  cout<<v.capacity()<<"\n"<<v.size()<<"\n"<<v[0];
  return 0;
}

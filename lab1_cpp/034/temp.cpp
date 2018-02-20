#include<list>
#include<iostream>
#include<vector>

using namespace std;

int main(){

  vector<vector<vector<int> > > v_list;
  vector<vector<int> > v1_list;
  vector<int> v2_list;

  v2_list.push_back(2);
  v1_list.push_back(2);
  v_list.push_back(v1_list);

  
  cout<<v_list.size()<<v_list.capacity();  
  //  v_list.push_back({0,1});

  cout<<v_list[0][0][0];
  return 0;
}


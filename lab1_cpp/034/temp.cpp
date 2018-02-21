#include<list>
#include<iostream>
#include<vector>

using namespace std;

int main(){

  int arr[]={16, 2, 3, 5};
  vector<int> vec(arr, arr+4);

for (int i=0;i<vec.size();i++){
  cout<<vec[i]<<"  ";
  }


  return 0;
}


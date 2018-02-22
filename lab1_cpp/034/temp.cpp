#include<list>
#include<iostream>
#include<vector>

using namespace std;

int main(){

  int arr[]={16, 2, 3, 5};
  vector<int> vec(arr, arr+4);

  int arr2[]={161, 22, 3, 5};
  vector<int> vec2(arr2, arr2+4);

for (int i=0;i<vec.size();i++){
  cout<<vec[i]<<"  ";
  }

for (int i=0;i<vec2.size();i++){
  cout<<vec2[i]<<"  ";
  }


 vector<int> vector(vec.begin()+2,vec.begin()+4);

 /*
 if (vector(vec.begin()+2,vec.begin()+4) == vector(vec2.begin()+2,vec2.begin()+4){
   cout<<"same";
   }
 */
  return 0;
}


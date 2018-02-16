#include<list>
#include<iostream>

using namespace std;

void f(int &temp_count){
  temp_count=temp_count+4;
}
int main(){

  list<int> listOfNumbers;

  //Inserting elements at end in list
  listOfNumbers.push_back(5);
  listOfNumbers.push_back(6);

  //Inserting elements at front in list
  listOfNumbers.push_front(2);
  listOfNumbers.push_front(1);

  for(list<int>::iterator list_iter = listOfNumbers.begin();
      list_iter != listOfNumbers.end(); list_iter++)
    {
      cout<<*list_iter<<endl;
    }


  int temp_count=0;
  f(temp_count);
  cout<<temp_count;

  return 0;
}


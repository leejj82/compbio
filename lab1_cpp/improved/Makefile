all : lab01_01 lab01_02 lab01_03

lab01_01 : lab01_01.o
	g++ -o lab01_01 lab01_01.o 

lab01_01.o : lab01_01.cpp
	g++ -c lab01_01.cpp

lab01_02 : lab01_02.o
	g++ -o lab01_02 lab01_02.o 

lab01_02.o : lab01_02.cpp
	g++ -c lab01_02.cpp

lab01_03 : lab01_03.o
	g++ -o lab01_03 lab01_03.o 

lab01_03.o : lab01_03.cpp
	g++ -c lab01_03.cpp

temp : temp.o
	g++ -o temp temp.o 

temp.o : temp.cpp
	g++ -c temp.cpp





clean :
	rm lab01_01 lab01_02 lab01_03 *.o temp 


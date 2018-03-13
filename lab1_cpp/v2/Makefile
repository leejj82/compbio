all : lab01_all

lab01_all : lab01_all.o
	g++ -o lab01_all lab01_all.o 

lab01_all.o : lab01_all.cpp
	g++ -c lab01_all.cpp

temp : temp.o
	g++ -o temp temp.o 

temp.o : temp.cpp
	g++ -c temp.cpp

clean :
	rm *.o temp lab01_all


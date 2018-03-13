all : lab01_01 lab01_01_sample lab01_01_debug

lab01_01 : lab01_01.o
	g++ -o lab01_01 lab01_01.o 

lab01_01_sample : lab01_01_sample.o
	g++ -o lab01_01_sample lab01_01_sample.o 

lab01_01_debug : lab01_01_debug.o
	g++ -o lab01_01_debug lab01_01_debug.o 

lab01_01.o : lab01_01.cpp
	g++ -c lab01_01.cpp

lab01_01_sample.o : lab01_01.cpp
	g++ -DSAMPLE -c lab01_01.cpp -o lab01_01_sample.o

lab01_01_debug.o : lab01_01.cpp
	g++ -DDEBUG -c lab01_01.cpp -o lab01_01_debug.o

clean :
	rm lab01_01 lab01_01_sample lab01_01_debug lab01_01*.o


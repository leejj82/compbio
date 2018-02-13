all : lab01 lab01_sample lab01_debug

lab01 : lab01.o
	g++ -o lab01 lab01.o 

lab01_sample : lab01_sample.o
	g++ -o lab01_sample lab01_sample.o 

lab01_debug : lab01_debug.o
	g++ -o lab01_debug lab01_debug.o 

lab01.o : lab01.cpp
	g++ -c lab01.cpp

lab01_sample.o : lab01.cpp
	g++ -DSAMPLE -c lab01.cpp -o lab01_sample.o

lab01_debug.o : lab01.cpp
	g++ -DDEBUG -c lab01.cpp -o lab01_debug.o

clean :
	rm lab01 lab01_sample lab01_debug lab01*.o


all : lab1 lab1_sample lab1_debug

lab1 : lab1.o
	g++ -o lab1 lab1.o 

lab1_sample : lab1_sample.o
	g++ -o lab1_sample lab1_sample.o 

lab1_debug : lab1_debug.o
	g++ -o lab1_debug lab1_debug.o 

lab1.o : lab1.cpp
	g++ -c lab1.cpp

lab1_sample.o : lab1.cpp
	g++ -DSAMPLE -c lab1.cpp -o lab1_sample.o

lab1_debug.o : lab1.cpp
	g++ -DDEBUG -c lab1.cpp -o lab1_debug.o

clean :
	rm lab1 lab1_sample lab1_debug lab1*.o


all : lab01_02 lab01_02_sample lab01_02_debug

lab01_02 : lab01_02.o
	g++ -o $@ lab01_02.o 

lab01_02_sample : lab01_02_sample.o
	g++ -o $@ lab01_02_sample.o 

lab01_02_debug : lab01_02_debug.o
	g++ -o $@ lab01_02_debug.o 

lab01_02.o : lab01_02.cpp
	g++ -c lab01_02.cpp

lab01_02_sample.o : lab01_02.cpp
	g++ -DSAMPLE -c lab01_02.cpp -o $@

lab01_02_debug.o : lab01_02.cpp
	g++ -DDEBUG -c lab01_02.cpp -o $@

temp : temp.o
	g++ -o temp temp.o

temp.o: temp.cpp
	g++ -c temp.cpp

clean :
	rm lab01_02 lab01_02_sample lab01_02_debug lab01_02*.o


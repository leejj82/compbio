all : lab01_03 lab01_03_sample lab01_03_debug

lab01_03 : lab01_03.o
	g++ -o $@ lab01_03.o 

lab01_03_sample : lab01_03_sample.o
	g++ -o $@ lab01_03_sample.o 

lab01_03_debug : lab01_03_debug.o
	g++ -o $@ lab01_03_debug.o 

lab01_03.o : lab01_03.cpp
	g++ -c lab01_03.cpp

lab01_03_sample.o : lab01_03.cpp
	g++ -DSAMPLE -c lab01_03.cpp -o $@

lab01_03_debug.o : lab01_03.cpp
	g++ -DDEBUG -c lab01_03.cpp -o $@

temp : temp.o
	g++ -o temp temp.o

temp.o: temp.cpp
	g++ -c temp.cpp

clean :
	rm lab01_03 lab01_03_sample lab01_03_debug lab01_03*.o


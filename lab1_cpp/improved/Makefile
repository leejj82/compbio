all : lab01_improved

lab01_improved : lab01_improved.o
	g++ -o lab01_improved lab01_improved.o 

lab01_improved.o : lab01_improved.cpp
	g++ -c lab01_improved.cpp

lab01_01_improved : lab01_01_improved.o
	g++ -o lab01_01_improved lab01_01_improved.o 

lab01_01_improved.o : lab01_01_improved.cpp
	g++ -c lab01_01_improved.cpp
temp : temp.o
	g++ -o temp temp.o 

temp.o : temp.cpp
	g++ -c temp.cpp

clean :
	rm *.o temp 


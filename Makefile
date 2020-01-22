Alignment: Alignment.o main.cpp
	g++ -o Alignment -g main.cpp Alignment.o

Alignment.o: Alignment.h Alignment.cpp
	g++ -c -g Alignment.cpp

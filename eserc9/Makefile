CC = g++
CFLAGS = -Wall -O3 --std=c++11

es9.exe : es9.o random.o eserc9classes.h eserc9functions.h
	$(CC) random.o es9.o -o es9.exe
es9.o : es9.cpp
	$(CC) -c es9.cpp -o es9.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es9.exe

CC = g++
CFLAGS = -Wall -O3 --std=c++11

es2.exe : es2.o random.o
	$(CC) random.o es2.o -o es2.exe
es2.o : es2.cpp
	$(CC) -c es2.cpp -o es2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o es2.exe data.out
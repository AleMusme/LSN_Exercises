CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1.exe : es1.o random.o eserc1classes.h eserc1functions.h
	$(CC) random.o es1.o -o es1.exe
es1.o : es1.cpp
	$(CC) -c es1.cpp -o es1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es1.exe

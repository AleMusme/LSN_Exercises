CC = g++
CFLAGS = -Wall -O3 --std=c++11

es3.exe : es3.o random.o 
	$(CC) random.o es3.o -o es3.exe
es3.o : es3.cpp 
	$(CC) -c es3.cpp -o es3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o es3.exe 
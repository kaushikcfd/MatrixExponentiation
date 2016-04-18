all: main

main: main.o
		g++ -g -std=c++11 main.o -o test -llapacke -lblas

main.o: main.cpp
		g++ -g -std=c++11 -c main.cpp -llapacke -lblas

clean:
		rm -rf *.o test *.gnu *.temp *.jpg load *.png *.dat *~


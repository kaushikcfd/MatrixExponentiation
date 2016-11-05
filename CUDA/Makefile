all: main

main: main.o
	nvcc main.o -o test -pg
main.o: main.cu
	nvcc -c main.cu -pg
clean:
	rm -rf *.o test *.dat

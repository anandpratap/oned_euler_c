all:
	icpc -openmp -O3 -Wall -wd588 -o main main.cpp 
tune:
	icpc -openmp -O3 -debug -Wall -vec-report1 -wd588 -g -o main main.cpp

clean:
	rm main
mic:
	icpc -openmp -O3 -mmic -Wall -vec-report1 -wd588 -g -o main main.cpp
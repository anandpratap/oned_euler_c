all:
	icpc -openmp -O3 -Wall -wd588 io.cpp base.cpp main.cpp -lm -o main 
tune:
	icpc -openmp -O3 -debug -Wall -vec-report1 -wd588 -g io.cpp base.cpp main.cpp -lm -o main  

clean:
	rm main
mic:
	icpc -openmp -O3 -mmic -Wall -wd588 io.cpp base.cpp main.cpp -lm -o main 
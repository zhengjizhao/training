all: test.x 

test.x: test.cpp 
	CC -dynamic -xMIC-AVX512 test.cpp -o test.x  -qopenmp -qopt-report=5 -qopt-report-phase=openmp,vec,loop -qopt-report-file=test.optrpt -std=c++17 -O3 -g -debug inline-debug-info -Wl,-rpath=/opt/gcc/7.2.0/snos/lib64

clean: 
	@rm -fr *.o *.x *.x.*


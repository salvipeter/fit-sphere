all: test

DEFINES= #-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
INCLUDES=-I/usr/include/eigen3
LIBS= #-lblas -llapacke
CXXFLAGS=-Wall -pedantic -std=c++17 -O3 $(INCLUDES) $(DEFINES)

test: test.o fit-sphere.o
	g++ -o $@ $^ $(LIBS)

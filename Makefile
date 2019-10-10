all: test

CXXFLAGS=-Wall -pedantic -std=c++17 -I/usr/include/eigen3

test: test.o fit-sphere.o
	g++ -o $@ $^

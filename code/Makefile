CC = gcc
CXX = g++
LIBS = -lrt
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 -mcmodel=large
TARGET = ./all

all: index CSP2Hop QHL
index: index.cpp
	$(CXX) $(CFLAGS) -o index index.cpp
	
CSP2Hop: CSP2Hop.cpp
	$(CXX) $(CFLAGS) -o CSP2Hop CSP2Hop.cpp
	
QHL: QHL.cpp
	$(CXX) $(CFLAGS) -o QHL QHL.cpp

.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm SOR shortestPath SRH Statistic 2> /dev/null || true

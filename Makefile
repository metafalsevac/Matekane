CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
HEADERS  := $(wildcard src/*.hpp)

all: matekane

matekane: Matekane.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ Matekane.cpp

tests/run_tests: tests/test_main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ tests/test_main.cpp

test: tests/run_tests
	./tests/run_tests

clean:
	rm -f matekane tests/run_tests

.PHONY: all test clean

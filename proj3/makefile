# compiler
CC=g++
# compiler flags
CFLAGS=-O3 -Wall -std=c++17 -fopenmp
# linker flags
LDFLAGS=
# source files
SOURCES=tsp-mpi.cpp
# object files
OBJECTS=$(SOURCES:.cpp=.o)
# target executable
EXECUTABLE=tsp

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SOURCES) -o $@


clean:
	rm -rf *o $(EXECUTABLE)

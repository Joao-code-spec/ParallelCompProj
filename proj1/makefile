# compiler
CC=g++
# compiler flags
CFLAGS=-c -Wall
# linker flags
LDFLAGS=
# source files
SOURCES=tsp.cpp
# object files
OBJECTS=$(SOURCES:.cpp=.o)
# target executable
EXECUTABLE=tsp

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o $(EXECUTABLE)

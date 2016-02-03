# Objects
TARGET = test
SOURCES = test.cpp Generator.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Compiler
CC = g++

# Flags
CFLAGS = -c -g
LFLAGS = -g

########################################

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o $@ 

%.o : %.cpp *.h
	$(CC) $< $(CFLAGS)

clean:
	-rm -f *.o *~ $(TARGET)

CC=gcc
CFLAGS=-c -Wall -O3 -march=native
LDFLAGS=-lm -lpthread -O3 -march=native
SOURCES=name.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=name

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm $(OBJECTS)
	rm $(EXECUTABLE)

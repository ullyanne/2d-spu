# Macros:
CPP = g++
CFLAGS = -g -O3 -w
OBJECTS = brkgarls.o Packing.o Decoder.o random_manager.o

# Targets:
all: brkgarls

brkgarls: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o brkgarls -lm

Decoder.o: Decoder.cpp Decoder.h
	$(CPP) $(CFLAGS) -c Decoder.cpp

brkgarls.o: brkgarls.cpp
	$(CPP) $(CFLAGS) -c brkgarls.cpp

Packing.o: Packing.cpp Packing.h
	$(CPP) $(CFLAGS) -c -g Packing.cpp

random_manager.o : random_manager.cpp random_manager.h
	$(CPP) $(CFLAGS) -c random_manager.cpp

# Remove:
clean:
	rm -f brkgarls $(OBJECTS)

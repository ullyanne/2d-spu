# Macros:
CPP = g++
CFLAGS = -g -O0 -fopenmp -Wall
OBJECTS = SampleDecoder.o samplecode.o Packing.o VirtualLayersDecoder.o

# Targets:
all: samplecode

samplecode: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o samplecode -lm

SampleDecoder.o: SampleDecoder.cpp SampleDecoder.h
	$(CPP) $(CFLAGS) -c SampleDecoder.cpp

VirtualLayersDecoder.o: VirtualLayersDecoder.cpp VirtualLayersDecoder.h
	$(CPP) $(CFLAGS) -c VirtualLayersDecoder.cpp

samplecode.o: samplecode.cpp
	$(CPP) $(CFLAGS) -c samplecode.cpp

Packing.o: Packing.cpp Packing.h
	$(CPP) $(CFLAGS) -c -g Packing.cpp

# Remove:
clean:
	rm -f samplecode $(OBJECTS)

# Macros:
CPP = g++
CFLAGS = -g -O3 -fopenmp -w
OBJECTS = samplecode.o Packing.o VirtualLayersDecoder.o random_manager.o

# Targets:
all: samplecode

samplecode: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o samplecode -lm

# SampleDecoder.o: SampleDecoder.cpp SampleDecoder.h
# 	$(CPP) $(CFLAGS) -c SampleDecoder.cpp

VirtualLayersDecoder.o: VirtualLayersDecoder.cpp VirtualLayersDecoder.h
	$(CPP) $(CFLAGS) -c VirtualLayersDecoder.cpp

samplecode.o: samplecode.cpp
	$(CPP) $(CFLAGS) -c samplecode.cpp

Packing.o: Packing.cpp Packing.h
	$(CPP) $(CFLAGS) -c -g Packing.cpp

random_manager.o : random_manager.cpp random_manager.h
	$(CPP) $(CFLAGS) -c random_manager.cpp

# Remove:
clean:
	rm -f samplecode $(OBJECTS)

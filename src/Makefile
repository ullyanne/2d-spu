# Macros:
CPP = g++
CFLAGS = -O3 -fopenmp -Wall
OBJECTS = SampleDecoder.o samplecode.o 

# Targets:
all: samplecode

samplecode: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o samplecode -lm
        
SampleDecoder.o: SampleDecoder.cpp SampleDecoder.h
	$(CPP) $(CFLAGS) -c SampleDecoder.cpp

samplecode.o: samplecode.cpp 
	$(CPP) $(CFLAGS) -c samplecode.cpp

# Remove:
clean:
	rm -f samplecode $(OBJECTS)

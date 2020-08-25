CC=g++ -std=c++11
CFLAGS =-c -g -Wall -DHAVE_INLINE -O0 -fexceptions -fopenmp
IFLAGS = -I /usr/include -I include -I /usr/local/include
LDFLAGS := -Xlinker -rpath -Xlinker /usr/local/lib -L /usr/local/lib -lgsl -lgslcblas -lCCfits -lcfitsio -fopenmp

#Common source files
CSOURCES= Img.cpp Model.cpp Sampler.cpp Writer.cpp GslRandom.cpp

#executable definitions
SOURCES=main.cpp $(CSOURCES)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)  

%.o: %.cpp
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *~ core */*.o */*~ *.op */*.op $(EXECUTABLE)

CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs   ibex)

CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-unused-variable -fopenmp -frounding-math -std=c++11 -I/usr/local/include/eigen -I/usr/local/libigl/include/igl

EXEC=test_dabbene test_kharitonov

all: $(EXEC)

test_kharitonov: kharitonov.o test_kharitonov.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)
	
test_dabbene: dabbene.o test_dabbene.o
	g++ -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	g++ -o $@ -c $< $(CXXFLAGS) $(LIBS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
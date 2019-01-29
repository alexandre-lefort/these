SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs  ibex)


#feq ($(DEBUG), yes)
#CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall -frounding-math -ffloat-store -fopenmp -Wno-unknown-pragmas -Wno-unused-variable
#CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -pedantic -frounding-math -ffloat-store  -Waddress -Warray-bounds  -Wc++11-compat  -Wchar-subscripts   -Wenum-compare    -Wcomment  -Wformat    -Wmain -Wmaybe-uninitialized  -Wmissing-braces -Wnonnull     -Wparentheses    -Wreorder  -Wreturn-type   -Wsequence-point  -Wsign-compare  -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunused-function -Wunused-label    -Wunused-value    -Wunused-variable -Wvolatile-register-var 
#else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-unused-variable -fopenmp -frounding-math -std=c++11
#endif

all: $(BINS)

% :	%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(BINS)

CC   = gcc-
CXX  = g++

CFLAGS = -std=c++11 -Wall -lfftw3 -lm

NVCCFLAGS = -std=c++11
LD_FLAGS =-std=c++11

EXEC = HMC

SRC_CPP = Main.cpp Monte_carlo.cpp stastics.cpp

OBJS_CPP = $(SRC_CPP:.cpp=.o)


all: $(EXEC)

$(EXEC): $(OBJS_CPP)
	$(CXX) -o $(EXEC) $(LD_FLAGS) $(OBJS_CPP)

clean:
	-rm *.o
	-rm *~
	-rm Make.log

distclean: clean
	-rm $(EXEC)

.SUFFIXES: .o .cpp

.cpp.o:
	$(CXX) $(CFLAGS) $(DEFS) $(INCS) -c $<

CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-Wall -std=c++17 -I/usr/include/eigen3 -march=native -fno-math-errno
LDFLAGS=-lm -lz


SRCS=cestimator.cpp support.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: cestimator

cestimator: $(OBJS)
	$(CXX) $(LDFLAGS) -o cestimator $(OBJS) $(LDLIBS) 

cestimator.o: cestimator.cpp support.hpp

support.o: support.hpp support.cpp

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) cestimator


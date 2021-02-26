# Load the common configuration file
#include ../../config.mk
cxx=g++-mp-10
cflags=-Wall -ansi -pedantic -O3 -fopenmp

objs=common.o ae_puzzle.o letter.o adj_table.o
src=$(patsubst %.o,%.cc,$(objs))
execs=ae_test

all: $(execs)
	
executables: $(execs)

depend: $(src)
	$(cxx) -MM $(src) >Makefile.dep

-include Makefile.dep

ae_test: ae_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all depend

include ../../config.mk

lflags=-llapack -lblas

# Lists of files to be built
objs=strassen.o standard.o extra.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mat_mul_test mat_mul_time

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

mat_mul_test: mat_mul_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lflags)

mat_mul_time: mat_mul_time.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lflags)

.PHONY: clean depend

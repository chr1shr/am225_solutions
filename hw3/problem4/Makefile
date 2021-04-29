# Include the configuration file
include ../../config.mk

# Use C++11 language standard
cflags+=-std=c++11

# Lists of files to be built
objs=helmholtz.o
src=$(patsubst %.o,%.cc,$(objs))
execs=helmholtz_test helmholtz_time

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(fftw_iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $<

helmholtz_test: helmholtz_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(lp_lflags) $(fftw_lflags)

helmholtz_time: helmholtz_time.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(lp_lflags) $(fftw_lflags)

.PHONY: clean depend

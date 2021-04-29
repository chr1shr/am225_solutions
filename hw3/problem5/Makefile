# Include the configuration file
include ../../config.mk

cflags+=-std=c++11

# Lists of files to be built
objs=poisson_fft.o poisson_schur.o subdomain.o conj_grad.o
src=$(patsubst %.o,%.cc,$(objs))
execs=schur_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(cflags) $(fftw_iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(fftw_iflags) -c $< $(fftw_lflags)

schur_test: schur_test.cc $(objs)
	$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(lflags) $(fftw_lflags) $(lp_lflags)

.PHONY: clean depend

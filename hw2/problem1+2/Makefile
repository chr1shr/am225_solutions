# Include the configuration file
include ../../config.mk

# Lists of files to be built
heads=brusselator.hh osc.hh
objs=sol_adapt.o sol_fsal.o sol_ckr.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=ode_test ode_conv

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

ode_test: ode_test.cc $(objs) $(heads)
	$(cxx) $(cflags) -o $@ $^

ode_conv: ode_conv.cc $(objs) $(heads)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend

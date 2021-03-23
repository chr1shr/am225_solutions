# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=sol_base.o sol_sun.o sol_rk4.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=ode_conv ode_test ga_test

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

ode_conv: ode_conv.cc $(objs) brusselator.hh
	$(cxx) $(cflags) -o $@ $^

ode_test: ode_test.cc $(objs) brusselator.hh
	$(cxx) $(cflags) -o $@ $^

ga_test: ga_test.cc $(objs) galaxy.hh
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend

# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=sol_adapt.o sol_fsal.o sol_ckr.o swarm.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=swarm_test

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

swarm_test: swarm_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend

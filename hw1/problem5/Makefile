# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=diffuse.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=diff_speed diff_stdev

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

diff_speed: diff_speed.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

diff_stdev: diff_stdev.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend

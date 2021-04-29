# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=lax_wend.o
src=$(patsubst %.o,%.cc,$(objs))
execs=lw_solve lw_conv

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

lw_solve: lw_solve.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

lw_conv: lw_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend

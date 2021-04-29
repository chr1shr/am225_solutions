# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=squircle.o conj_grad.o quadrat.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=sq_test sq_conv

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
	$(cxx) $(cflags) -c $<

sq_test: sq_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

sq_conv: sq_conv.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend

# Include the configuration file
include ../../config.mk

adapt_dir=../problem5

iflags=-I$(adapt_dir)

# Lists of files to be built
objs=sol_rk4.o sol_base.o
extra_objs=$(adapt_dir)/sol_ckr.o $(adapt_dir)/sol_fsal.o $(adapt_dir)/sol_adapt.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=dc_conv dc_test

all:
	$(MAKE) -C $(adapt_dir)
	$(MAKE) $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

dc_test: dc_test.cc discont.hh $(objs) $(extra_objs)
	$(cxx) $(cflags) $(iflags) -o $@ $< $(objs) $(extra_objs)

dc_conv: dc_conv.cc discont.hh $(objs) $(extra_objs)
	$(cxx) $(cflags) $(iflags) -o $@ $< $(objs) $(extra_objs)

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

.PHONY: clean depend

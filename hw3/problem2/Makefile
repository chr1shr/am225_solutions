# Path to main RBF code
rbfp=../../../lectures/unit2/lec07_example

include ../../config.mk

# Lists of files to be built
objs=rbf_hilbert.o rbf.o conj_grad.o lp_solve.o file_output.o
src=$(patsubst %.o,%.cc,$(objs))
execs=rbfh_test rbfh_time

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
	$(cxx) $(cflags) $(iflags) -c $<

rbfh_test: rbfh_test.cc $(objs)
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

rbfh_time: rbfh_time.cc $(objs)
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend

include ../../config.mk

execs=tree_test tree_count bh_trees
objs=tree.o
src=$(patsubst %.o,%.cc,$(objs))

all: $(execs)
-include Makefile.dep

depend:
	$(cxx) $(iflags) -MM $(src) > Makefile.dep

clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

tree_test: tree_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

tree_count: tree_count.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

bh_trees: bh_trees.cc
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend

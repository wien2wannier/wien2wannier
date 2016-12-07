### wien2wannier/test/Makefile
###
###    wien2wannier/test/Makefile
###
### Copyright 2016 Elias Assmann

-include ../make.sys

SUBDIRS := PoC

.PHONY: all clean $(SUBDIRS) dist

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		rm -f $$dir/Makefile.orig; \
	done

distclean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done


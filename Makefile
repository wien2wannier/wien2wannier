### wien2wannier/Makefile
###
###    wien2wannier main Makefile
###
### Copyright 2013-2015 Elias Assmann

VERSION := $(shell git describe)
ifeq "$(VERSION)" ""
VERSION = $(lastword '$version: v1.0.0-125-g9f7266f$')
endif

SIMPLE      := SRC_trig doc
REALCOMPLEX := SRC_w2w SRC_wplot

SUBDIRS := $(SIMPLE) $(REALCOMPLEX)

.PHONY: all clean $(SUBDIRS) dist

all: $(SUBDIRS)

$(REALCOMPLEX):
	$(MAKE) -C $@ real
	$(MAKE) -C $@ complex

$(SIMPLE):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
		rm -f $$dir/Makefile.orig; \
	done

distclean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir distclean; \
	done

%/Makefile.orig: %/Makefile
	perl -pe 's/^\#.orig\#//' $^ >$@

Morig := $(addsuffix /Makefile.orig,SRC_w2w SRC_wplot SRC_trig)

## Make a tarball for extracting directly in Wien2k root folder and
## bundle with INSTALL, compile_wien2wannier.sh
##
## This should perhaps strip the time-stamps.
wien-dist: dir     = wien2wannier-$(VERSION)-wiendist
wien-dist: tarname = wien2wannier-$(VERSION)-expand-in-wienroot
wien-dist: scripts = $(notdir $(wildcard SRC/*))
wien-dist: w2wlinks= ../doc/wien2wannier_userguide.pdf ../COPYING \
	   	     ../README ../NEWS ../doc/CHEATSHEET
wien-dist: distclean doc/wien2wannier_userguide.pdf $(Morig)
	mkdir $(dir); \
	cd $(dir); \
	ln -s -t. ../SRC* \
		../compile_wien2wannier.sh ../INSTALL ../WIEN-VERSION; \
	cp -t. SRC/*; \
	ln -s -tSRC_w2w $(w2wlinks); \
	tar --exclude-vcs -chf $(tarname).tar SRC* $(scripts); \
	rm SRC* $(scripts)

	tar --exclude-vcs -chzf $(dir).tar.gz $(dir)

	for l in $(notdir $(w2wlinks)); do rm SRC_w2w/$$l; done
	rm -rf $(dir) $(Morig)


## Time-stamp: <2016-02-09 14:14:51 assman@faepop36.tu-graz.ac.at>

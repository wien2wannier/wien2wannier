### wien2wannier/Makefile
###
###    wien2wannier main Makefile
###
### Copyright 2013-2016 Elias Assmann

-include make.sys

SHELL=/bin/bash

ifndef VERSION
   VERSION = $(lastword '$version: v1.0.0-209-ga0388eb$')
endif

SIMPLE      := SRC_trig doc test
REALCOMPLEX := SRC_w2w SRC_wplot

SUBDIRS := $(SIMPLE) $(REALCOMPLEX)

.PHONY: all clean distclean $(SUBDIRS) dist-tmp wien-tar wien-dist install

all: SRC_w2w SRC_wplot SRC_trig

$(REALCOMPLEX):
	$(MAKE) -C $@ real
	$(MAKE) -C $@ complex

$(SIMPLE):
	$(MAKE) -C $@

test: target-dir = $(WIENROOT_TEST)
test: install

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


### ‘install’ target
###
### This is a bare-bones install procedure for updating wien2wannier
### in an existing Wien2k directory.  Use with caution.
target-dir ?= $(WIENROOT)
install: files = $(shell find SRC_{w2w,wplot,trig} -type f -perm /a+x)
install: all
	install -t$(target-dir) $(files)


### Distribution targets

## dist variables
Morig    := $(addsuffix /Makefile.orig,SRC_w2w SRC_wplot SRC_trig)
dist-dir := wien2wannier-$(VERSION)-wiendist
tarname  := wien2wannier-$(VERSION)-expand-in-wienroot
scripts  := $(notdir $(wildcard SRC/*))
w2wlinks := ../doc/wien2wannier_userguide.pdf \
	    ../COPYING ../README ../NEWS ../doc/CHEATSHEET

## Make a tarball for extracting directly in Wien2k root folder
##
## This should perhaps strip the time-stamps.
dist-tmp: distclean doc/wien2wannier_userguide.pdf $(Morig)
	mkdir $(dist-dir); \
	cd $(dist-dir); \
	ln -s -t. ../SRC* \
		../compile_wien2wannier.sh ../INSTALL ../WIEN-VERSION; \
	cp -t. SRC/*; \
	ln -s -tSRC_w2w $(w2wlinks); \
	tar --exclude-vcs -chf $(tarname).tar SRC* $(scripts);

## Make tarball and clean up
wien-tar: dist-tmp
	mv $(dist-dir)/$(tarname).tar .

	cd $(dist-dir); \
	rm SRC* $(scripts)

	for l in $(notdir $(w2wlinks)); do rm SRC_w2w/$$l; done
	rm -rf $(dist-dir) $(Morig)

## Make tarball and bundle with INSTALL, compile_wien2wannier.sh
wien-dist: dist-tmp
	cd $(dist-dir); \
	rm SRC* $(scripts)

	tar --exclude-vcs -chzf $(dist-dir).tar.gz $(dist-dir)

	for l in $(notdir $(w2wlinks)); do rm SRC_w2w/$$l; done
	rm -rf $(dist-dir) $(Morig)


## Time-stamp: <2016-07-26 16:15:18 assman@faepop71.tu-graz.ac.at>

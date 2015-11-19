### wien2wannier/Makefile
###
###    wien2wannier main Makefile
###
### Copyright 2013-2015 Elias Assmann
###
### $Id: Makefile 422 2015-07-01 08:24:37Z assmann $

svn-rev := r$(lastword '$Rev: 422 $')

VERSION := $(svn-rev)

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

dist: dir     = wien2wannier-$(VERSION)
dist: scripts = $(notdir $(wildcard SRC/*))
dist: distclean doc/wien2wannier_userguide.pdf
	mkdir $(dir); \
	cd $(dir); \
	ln -s -t . ../SRC* ../doc/ ../COPYING \
           ../README ../NEWS ../Makefile; \
	ln -s ../make.sys make.sys.example; \
	cp ../WIEN-VERSION .

	tar --exclude-vcs -chzf $(dir).tar.gz $(dir)
	rm -rf $(dir) $(Morig)

old-dist: dir     = wien2wannier-$(VERSION)
old-dist: scripts = $(notdir $(wildcard SRC/*))
old-dist: distclean doc/wien2wannier_userguide.pdf $(Morig)
	mkdir $(dir); \
	cd $(dir); \
	ln -s -t . ../SRC* ../doc ../COPYING \
	   ../README ../INSTALL ../NEWS ../Makefile \
	   ../compile_wien2wannier.sh; \
	ln -s ../make.sys make.sys.example; \
	ln -s -t . SRC/*; \
	cp ../WIEN-VERSION .; \
	tar --exclude-vcs -chf wien2wannier.tar SRC* $(scripts); \
	# for f in $(scripts); do ln -s $$f `echo $$f | sed s/_lapw//`; done; \
	# ln -s w2wpara w2wcpara; ln -s wplotpara wplotcpara; \
	# tar -rf wien2wannier.tar $(scripts:_lapw=) w2wcpara wplotcpara; \
	# rm -f SRC* $(scripts) $(scripts:_lapw=) w2wcpara wplotcpara

	# tar --exclude-vcs -chzf $(dir).tar.gz $(dir)
	# rm -rf $(dir) $(Morig)

## Make a .tar for extracting directly in Wien2k root folder and
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


## Time-stamp: <2015-07-01 10:23:12 assman@faepop23.tu-graz.ac.at>

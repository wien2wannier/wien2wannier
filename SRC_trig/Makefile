### wien2wannier/SRC_trig/Makefile

####### Configuration variables #######

SHELL = /bin/sh
.SUFFIXES:

FFLAGS += $(FGEN) $(FOPT)
#.orig#DESTDIR:= ../

#.orig### Variables inherited from WIEN
#.orig### these will be written by ‘siteconfig’
#.orig#FC     	:= _FC_
#.orig#CC     	:= _CC_
#.orig#FOPT   	:= _FOPT_
#.orig#LDFLAGS	:= _LDFLAGS_
#.orig#R_LIBS 	:= _R_LIBS_

#.orig## $(MODOPT) should be the option that tells the compiler where to put
#.orig## .mod files
#.orig##
#.orig## gfortran: -J
#.orig## ifort:    -module
#.orig#ifneq "$(findstring gfortran,$(FC))" ""
#.orig#MODOPT ?= -J
#.orig#else ifneq "$(findstring ifort,$(FC))" ""
#.orig#MODOPT ?= -module
#.orig#else ifndef MODOPT
#.orig#$(error unknown compiler ``$(FC)'', please specify variables manually)
#.orig#endif

SRCDIR  := ../SRC_w2w
LIBDIR  := lib
UTILDIR := $(SRCDIR)/$(LIBDIR)
UTILOBJ := $(UTILDIR)/util_w2w.o

## site-specific configuration may overwrite the above
-include ../make.sys

#######################################

EXE_UTIL := write_win_backend findbands shifteig write_inwplot	\
	    rephase w2waddsp join_vectorfiles vec2ascii		\
	    convham

#.orig#EXE_SIMPLE := rhomb_in5 hex2rhomb rhomb2hex xyz2struct convert_in2	\
#.orig#	      plane clminter struct2cif Pinver fmax add_columns	sub_columns

EXE := $(EXE_UTIL) $(EXE_SIMPLE)

#.orig#EXE += str2xcr.exe

UTIL_USERS := $(addsuffix .o,$(EXE_UTIL))

all: $(EXE)
#.orig#	cp $(EXE) $(DESTDIR)

clean:
	rm -f *.mod *.o

distclean: clean
	rm -f $(EXE)

.PHONY: all clean distclean

####

$(UTILOBJ): $(SRCDIR)/util_w2w.F $(SRCDIR)/Makefile
	make -C $(SRCDIR) $(LIBDIR)/util_w2w.o

$(EXE_UTIL): %: %.o $(UTILOBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(EXE_SIMPLE): %: %.f
	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(UTIL_USERS): $(UTILOBJ)
$(UTIL_USERS): FFLAGS += -I$(UTILDIR)

#.orig#str2xcr.exe: str2xcr.f auxil.o
#.orig#	$(FC) $(FFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $<


## Time-stamp: <2016-02-04 14:23:25 assman@faepop36.tu-graz.ac.at>

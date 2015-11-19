### wien2wannier/SRC_wplot/Makefile
###
###    make           ... generate executable for the real-valued version
###    make complex   ... generate executable for the complex-valued version
###    make all       ... generate both the above
###    make clean     ... delete unnecessary files
###
### $Id: Makefile 409 2015-06-09 07:54:07Z assmann $


# FC ........... compiler name
# FOPT ......... optimization flags
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# DESTDIR ...... dir. where the executable should go (without trailing '/'!)
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )

SHELL = /bin/sh
.SUFFIXES:

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

C_LIBS = $(R_LIBS)

DESTDIR    := .
R_EXECNAME := wplot
C_EXECNAME := wplotc

# directories for ‘.o’ and ‘.mod’ files that are: common to real and
# complex; real versions; complex versions
O_DIR := lib
R_DIR := libr
C_DIR := libc

###############################################################################

FFLAGS  = $(FGEN) $(FOPT) -I$(O_DIR)
R_EXEC  = $(DESTDIR)/$(R_EXECNAME)
C_EXEC  = $(DESTDIR)/$(C_EXECNAME)

SRCDIR  := ../SRC_w2w
LIBDIR  := lib
UTILDIR := $(SRCDIR)/$(LIBDIR)
UTILOBJ := $(UTILDIR)/util_w2w.o

## site-specific configuration may overwrite the above
-include ../make.sys

#..............................................................................
#
#  convenience rules
#
# by WIEN convention, ‘real’ is the default target, and ‘clean’ does
# not remove executables
real:    $(R_EXEC)
complex: $(C_EXEC)

all: real complex

clean:
	rm -f $(objmod)

distclean: clean
	rm -f $(R_EXEC) $(C_EXEC)

.PHONY: real complex all clean distclean

#..............................................................................
#
#  Object files common to both REAL and COMPLEX version
#  (i.e., these are the same for REAL and COMPLEX)
#
OBJS = moduls.o wplot.o latgen.o rotdef.o locdef.o orth.o grdgen.o	\
       spcgen.o findmt.o auggen.o augpw.o auglo.o wavsph.o wavint.o	\
       trans.o bessel.o gbass.o outwin.o rint13.o ylm.o sphbes.o	\
       dvbes1.o diracout.o dergl.o inouh.o inth.o

#..............................................................................
#
#  These object files have REAL and COMPLEX versions
#
#  The two versions go to different directories so we can keep both
#  versions around.
#
R_OBJS = main.o
C_OBJS = $(R_OBJS)

#..............................................................................
#
#  Dependency information for object files
#
MOD_USERS = auggen.o auglo.o augpw.o dergl.o diracout.o findmt.o	\
	    grdgen.o inouh.o locdef.o outwin.o rint13.o rotdef.o	\
	    wavint.o wavsph.o main.o

UTIL_USERS = auggen.o augpw.o auglo.o bessel.o dergl.o dvbes1.o		\
	     grdgen.o gbass.o inouh.o inth.o latgen.o locdef.o main.o	\
	     moduls.o orth.o outwin.o rotdef.o spcgen.o sphbes.o	\
	     trans.o wavint.o wplot.o ylm.o

#..............................................................................
#
#  “Real” targets
#
objs   = $(foreach o, $(OBJS)  , $(O_DIR)/$o)
r_objs = $(foreach o, $(R_OBJS), $(R_DIR)/$o)
c_objs = $(foreach o, $(C_OBJS), $(C_DIR)/$o)

objmod = $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
         $(foreach s, o mod, $d/*.$s))

$(R_EXEC): $(objs) $(r_objs) $(UTILOBJ) Makefile
	$(FC) -o $(R_EXEC) $(objs) $(r_objs) $(UTILOBJ) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(objs) $(c_objs) $(UTILOBJ) Makefile
	$(FC) -o $(C_EXEC) $(objs) $(c_objs) $(UTILOBJ)  $(LDFLAGS) $(C_LIBS)

#..............................................................................
#
#  Pattern rules for object files
#
#  Real and complex need to be kept separate!
#
$(O_DIR)/%.o: %.f Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(O_DIR)
$(O_DIR)/%.o: %.F Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(O_DIR)

$(R_DIR)/%.o: %.f Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(R_DIR)
$(R_DIR)/%.o: %.F Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(R_DIR)

$(C_DIR)/%.o: %.f Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(C_DIR) -D_COMPLEX_
$(C_DIR)/%.o: %.F Makefile
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(C_DIR) -D_COMPLEX_

#..............................................................................
#
#  Dependencies
#
$(UTILOBJ): $(SRCDIR)/util_w2w.F $(SRCDIR)/Makefile
	make -C $(SRCDIR) $(LIBDIR)/util_w2w.o

# Dependencies on ‘moduls.o’ and ‘util_w2w.o’; these files might have
# real/complex versions, or not
mod_users  = $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
		$(foreach o, $(MOD_USERS), $d/$o))
util_users = $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
		$(foreach o, $(UTIL_USERS), $d/$o))

$(mod_users): $(O_DIR)/moduls.o

$(util_users): $(UTILOBJ)
$(util_users): FFLAGS += -I$(UTILDIR)


## Time-stamp: <2015-06-08 15:56:16 assman@faepop23.tu-graz.ac.at>

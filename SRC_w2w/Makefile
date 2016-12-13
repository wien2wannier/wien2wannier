### wien2wannier/SRC_w2w/Makefile
###
###    make         ... generate executable for the real-valued version
###    make complex ... generate executable for the complex-valued version
###    make all     ... generate both the above
###    make clean   ... delete unnecessary files
###
### Copyright 2013-2016 Elias Assmann

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

DESTDIR	    = .
R_EXECNAME := w2w
C_EXECNAME := $(R_EXECNAME)c

# directories for ‘.o’ and ‘.mod’ files that are: common to real and
# complex; real versions; complex versions
O_DIR := lib
R_DIR := libr
C_DIR := libc

#############################################################################

FFLAGS+= $(FGEN) $(FOPT) -I$(O_DIR)
R_EXEC = $(DESTDIR)/$(R_EXECNAME)
C_EXEC = $(DESTDIR)/$(C_EXECNAME)

## site-specific configuration may overwrite the above
-include ../make.sys

#............................................................................
#
#  convenience rules
#
# by WIEN convention, ‘real’ is the default target, and ‘clean’ does
# not remove executables
real:    $(R_EXEC)
complex: $(C_EXEC)

all: real complex

clean:
	rm -f $(objs) $(r_objs) $(c_objs) $(O_DIR)/*.mod $(R_DIR)/*.mod $(C_DIR)/*.mod

distclean: clean
	rm -df $(R_EXEC) $(C_EXEC) f $(O_DIR) $(R_DIR) $(C_DIR)

.PHONY: real complex all clean distclean

#........................................................................
#
#  These object files have REAL and COMPLEX versions
#
#  The two versions go to different directories so we can keep both
#  versions around.
#
R_OBJS = modules_rc.o
C_OBJS = $(R_OBJS)

#........................................................................
#
#  Object files that depend on REAL/COMPLEX indirectly
#
#  We must exclude those objects that are already named in R_OBJS
#
RC_OBJS = $(filter-out $(R_OBJS), $(modrc_users))

#........................................................................
#
#  Object files common to both REAL and COMPLEX version
#  (i.e., these are the same for REAL and COMPLEX)
#
OBJS = $(filter-out $(R_OBJS) $(RC_OBJS), \
	$(patsubst %.F, %.o, $(wildcard *.F)) \
	$(patsubst %.f, %.o, $(wildcard *.f)))

#........................................................................
#
#  Dependency information for object files
#

## util_w2w.F:   utilities used in SRC_w2w, SRC_wplot, SRC_trig
util_users   = mod_w2w.o modules.o modules_rc.o main.o
## mod_w2w.f:    modules and procedures shared by SRC_w2w and SRC_wplot
modw2w_users = modules.o modules_rc.o main.o
## modules.f:    modules and procedures for SRC_w2w
mod_users    = modules_rc.o main.o
## modules_rc.f: modules and procedures for SRC_w2w -- real/complex
modrc_users  = main.o

#........................................................................
#
#  “Real” targets
#
objs   = $(foreach o, $(OBJS)             , $(O_DIR)/$o)
r_objs = $(foreach o, $(RC_OBJS) $(R_OBJS), $(R_DIR)/$o)
c_objs = $(foreach o, $(RC_OBJS) $(C_OBJS), $(C_DIR)/$o)

$(R_EXEC): $(objs) $(r_objs) Makefile
	$(FC) -o $(R_EXEC) $(objs) $(r_objs) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(objs) $(c_objs) Makefile
	$(FC) -o $(C_EXEC) $(objs) $(c_objs) $(LDFLAGS) $(C_LIBS)

#........................................................................
#
#  Pattern rules for object files
#
#  Real and complex need to be kept separate!
#
$(O_DIR)/%.o: %.[fF] Makefile
	@test -d $(O_DIR) || mkdir $(O_DIR)
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(O_DIR)

$(O_DIR)/%.o: %.c Makefile
	@test -d $(O_DIR) || mkdir $(O_DIR)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

$(R_DIR)/%.o: %.[fF] Makefile
	@test -d $(R_DIR) || mkdir $(R_DIR)
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(R_DIR)

$(C_DIR)/%.o: %.[fF] Makefile
	@test -d $(C_DIR) || mkdir $(C_DIR)
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(C_DIR) -D_COMPLEX_

#........................................................................
#
#  Dependencies
#

## Dependencies on ‘util_w2w.o’; these files might have real/complex
## versions, or not
$(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(util_users), $d/$o)): $(O_DIR)/util_w2w.o

## Dependencies on ‘mod_w2w.o’; these files might have real/complex
## versions, or not
$(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(modw2w_users), $d/$o)): $(O_DIR)/mod_w2w.o

## Dependencies on ‘modules.o’; these files might have real/complex
## versions, or not
$(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(mod_users), $d/$o)): $(O_DIR)/modules.o

## Dependencies on ‘modules_rc.o’; we know these files must have
## real/complex versions
$(foreach o, $(modrc_users), $(R_DIR)/$o): $(R_DIR)/modules_rc.o
$(foreach o, $(modrc_users), $(C_DIR)/$o): $(C_DIR)/modules_rc.o

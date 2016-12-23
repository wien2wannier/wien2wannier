### wien2wannier/SRC_wplot/Makefile
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

DESTDIR    := .
R_EXECNAME := wplot
C_EXECNAME := $(R_EXECNAME)c

# directories for ‘.o’ and ‘.mod’ files that are: common to real and
# complex; real versions; complex versions
O_DIR := lib
R_DIR := libr
C_DIR := libc

#########################################################################

FFLAGS += $(FGEN) $(FOPT) -I$(O_DIR)
R_EXEC  = $(DESTDIR)/$(R_EXECNAME)
C_EXEC  = $(DESTDIR)/$(C_EXECNAME)

## where to find the stuff we include from w2w
SRCDIR  := ../SRC_w2w
LIBDIR  := lib
UTILDIR := $(SRCDIR)/$(LIBDIR)
UTILOBJ := $(UTILDIR)/util_w2w.o
MODOBJ  := $(UTILDIR)/mod_w2w.o

## site-specific configuration may overwrite the above
-include ../make.sys

#........................................................................
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
#  Object files common to both REAL and COMPLEX version
#  (i.e., these are the same for REAL and COMPLEX)
#
OBJS = $(filter-out $(R_OBJS) $(RC_OBJS), \
	$(patsubst %.F, %.o, $(wildcard *.F)) \
	$(patsubst %.f, %.o, $(wildcard *.f)))

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
#  Dependency information for object files
#

## util_w2w.F:   utilities used in SRC_w2w, SRC_wplot, SRC_trig
util_users   = modules.o modules_rc.o wplot.o
## mod_w2w.f:    modules and procedures shared by SRC_w2w and SRC_wplot
modw2w_users = modules.o modules_rc.o modules_rc.o wplot.o
## modules.f:    modules and procedures for SRC_wplot
mod_users    = modules_rc.o wplot.o
## modules_rc.f: modules and procedures for SRC_wplot -- real/complex
modrc_users  = wplot.o

#........................................................................
#
#  “Real” targets
#
objs   = $(foreach o, $(OBJS)             , $(O_DIR)/$o)
r_objs = $(foreach o, $(RC_OBJS) $(R_OBJS), $(R_DIR)/$o)
c_objs = $(foreach o, $(RC_OBJS) $(C_OBJS), $(C_DIR)/$o)

$(R_EXEC): $(objs) $(r_objs) $(UTILOBJ) $(MODOBJ) Makefile
	$(FC) -o $(R_EXEC) $(objs) $(r_objs) $(UTILOBJ) \
		$(MODOBJ) $(LDFLAGS) $(R_LIBS)

$(C_EXEC): $(objs) $(c_objs) $(UTILOBJ) $(MODOBJ) Makefile
	$(FC) -o $(C_EXEC) $(objs) $(c_objs) $(UTILOBJ) \
		$(MODOBJ) $(LDFLAGS) $(C_LIBS)

#........................................................................
#
#  Pattern rules for object files
#
#  Real and complex need to be kept separate!
#
$(O_DIR)/%.o: %.[fF] Makefile
	@test -d $(O_DIR) || mkdir $(O_DIR)
	$(FC) $(FFLAGS) -c $< -o$@ $(MODOPT) $(O_DIR)

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
$(UTILOBJ): $(SRCDIR)/util_w2w.F $(SRCDIR)/Makefile
	make -C $(SRCDIR) $(LIBDIR)/util_w2w.o

$(MODOBJ): $(SRCDIR)/mod_w2w.f $(SRCDIR)/Makefile
	make -C $(SRCDIR) $(LIBDIR)/mod_w2w.o

# Dependencies on ‘modules.o’, ‘mod_w2w.o’, and ‘util_w2w.o’; these
# files might have real/complex versions, or not, so we prefix them
# with each of the repsective directories
mod_users    := $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(mod_users), $d/$o))

modw2w_users := $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(modw2w_users), $d/$o))

util_users   := $(foreach d, $(O_DIR) $(R_DIR) $(C_DIR), \
	$(foreach o, $(util_users), $d/$o))

$(mod_users):    $(O_DIR)/modules.o

$(modw2w_users): $(MODOBJ)
$(modw2w_users): FFLAGS += -I$(UTILDIR)

$(util_users):   $(UTILOBJ)
$(util_users):   FFLAGS += -I$(UTILDIR)

# Dependencies on ‘modules_rc.o’; we know these files must have
# real/complex versions
$(foreach o, $(modrc_users), $(R_DIR)/$o): $(R_DIR)/modules_rc.o
$(foreach o, $(modrc_users), $(C_DIR)/$o): $(C_DIR)/modules_rc.o

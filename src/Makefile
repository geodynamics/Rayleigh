#############################################################
#  Machine_Definitions contains the various machine-specific 
#  definitions for FFLAGS_*, LIB_*, and VAR_*
include Machine_Definitions

include object_list  #List of object files and associated variables

######################################################

MAKENP  = $(MAKE) --no-print-directory
.SUFFIXES: .o .F90 .f .F .c

# White-space-separated list of all acceptable goals for this makefile
VALIDGOALS = $(VERSIONS) rayleigh_exec all default clean clean_obj clean_exec .F90.o .cc.o

#Check to see if the  user has specified a valid goal
#GOALCHECK will be an empty string if  MAKECMDGOALS is not found within VALIDGOALS
GOALCHECK = $(findstring $(MAKECMDGOALS),$(VALIDGOALS))

ifeq ($(GOALCHECK), $(MAKECMDGOALS))
  #If the goal is valid, check to see if it happens to be
  #referring to a user-specified version
  VCHECK = $(findstring $(MAKECMDGOALS),$(VERSIONS))
  ifdef VCHECK
    # If so, remap the version variable to version name
    # Make will then proceed to execute the $(VERSION) 
    # recipe below.
    VERSION=$(VCHECK)
    VV := $(FFLAGS_$(VERSION))
  endif
else
  #If the goal is invalid, print out some explanatory information.
  $(info )
  $(info ----ERROR: invalid target specification----)
  $(info Specified Target            :   $(MAKECMDGOALS))
  $(info Valid Intrinsic Targets     :   all clean clean_obj clean_exec)
  $(info Valid User-defined Targets:   $(VERSIONS))
  $(info )
endif


default:
	@$(MAKENP) $(DEFAULT) VERSION=$(DEFAULT) PREFIX=$(PREFIX)

all: 
	@for i in $(VERSIONS) ; do \
		$(MAKENP) $$i VERSION=$$i PREFIX=$(PREFIX); \
		echo ' ' ; \
	done

rayleigh_exec: $(OBJ)
	$(FC) $(FFLAGS) -o  $(PROG) $(OBJ) $(LIB)

# the lines that contain
#    cp $(CUSTOMROOT)/*.F . || :
# are basically logical OR statements, but the ":" shell character is always True.
# So the end result is to suppress any and all errors from a failed cp command
$(VERSION):
	@echo ' ' 
	@echo '******************************************************************** ' 
	@echo '                 Now Building: '$(PREFIX)'/rayleigh'.$(VERSION) 
	@echo ' '
	@echo '                  '$(DESC_$(VERSION))	
	@echo ' '
	@$(MAKENP) clean_obj
	@cp ../Physics/Run_Parameters.F90 .
	$(eval git_hash := $(shell git log -1 --pretty=format:"%H --- %cd"))
	$(eval git_url := $(shell git config remote.origin.url))
	$(eval git_branch := $(shell git rev-parse --abbrev-ref HEAD))
	$(eval git_hash := $(shell git log -1 --pretty=format:"%H --- %cd"))
	$(eval build_date := $(shell date))
	$(eval build_machine := $(shell uname -a))
	$(eval build_dir := $(shell pwd))
	$(eval build_version := $(shell git describe --tags))
	$(eval compiler_version := $(shell $(FC) --version))
	$(eval compiler_location := $(shell which $(FC)))
	@sed -i'' "s|\@\@git hash\@\@|$(git_hash)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@git url\@\@|$(git_url)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@git branch\@\@|$(git_branch)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@date\@\@|$(build_date)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@machine\@\@|$(build_machine)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@directory\@\@|$(build_dir)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@custom dir\@\@|$(CUSTOMROOT)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@build version\@\@|$(VERSION)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@fflags\@\@|$(FFLAGS_$(VERSION))|g" Run_Parameters.F90
	@sed -i'' "s|\@\@library\@\@|$(LIB_$(VERSION))|g" Run_Parameters.F90
	@sed -i'' "s|\@\@FC version\@\@|$(compiler_version)|g" Run_Parameters.F90
	@sed -i'' "s|\@\@FC location\@\@|$(compiler_location)|g" Run_Parameters.F90
	@python adjust_Run_Parameters.py Run_Parameters.F90 Run_Parameters.$(VERSION).F90
	@cp Run_Parameters.$(VERSION).F90 Run_Parameters.F90
ifdef CUSTOMROOT
	@cp $(CUSTOMROOT)/*.F90 . 2> /dev/null || :
	@cp $(CUSTOMROOT)/*.F . 2> /dev/null || :
	@cp $(CUSTOMROOT)/*.c . 2> /dev/null || :
endif
	@$(MAKENP) rayleigh_exec FFLAGS="$(FFLAGS_$(VERSION))" PROG="compiled/rayleigh.$(VERSION)" LIB="$(LIB_$(VERSION))" PREFIX=$(PREFIX)

.F90.o :
	$(FC) $(FFLAGS) -c $<
.cc.o:
	$(CC) $(CFLAGS) -c *.cc

clean_exec:
	@echo "Cleaning up existing executables in src/build/compiled..."
	@rm -rf compiled/*
clean_obj:
	@echo "Cleaning up compiled object files..."
	@rm -f *.o *.mod
clean: 
	@$(MAKENP) clean_obj


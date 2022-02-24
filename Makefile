-include make.inc
SRC=src
BUILD=$(SRC)/build
INTERP=post_processing/interpolation

# LC_COLLATE is set because sort is locale dependent otherwise.
export LC_COLLATE=C

# make the CUSTOMROOT variable available to sub-make processes
export CUSTOMROOT

# running "make target=dbg" will only compile the specified target
target=all

rayleigh: prepare_directory
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) $(target)
	@echo ""
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Compilation is complete."
	@echo " "
	@echo " "
	@echo "Run 'make install' to install the rayleigh executables into: "
	@echo $(PREFIX)"/bin"
	@echo ""
fdeps: prepare_directory
	$(MAKE) --no-print-directory --directory=$(BUILD) fdeps
	cp $(BUILD)/Makefile.fdeps $(SRC)/.
.PHONY: prepare_directory
prepare_directory:
	@mkdir -p $(BUILD)/compiled
	@cp $(SRC)/Makefile.fdeps $(BUILD)/.
	@cp $(SRC)/Parallel_Framework/*.F90 $(BUILD)/.
	@cp $(SRC)/Data_Structures/*.F90 $(BUILD)/.
	@cp $(SRC)/Math_Layer/*.F90 $(BUILD)/.
	@cp $(SRC)/Plugins/*.F90 $(BUILD)/.  2>/dev/null || :
	@cp $(SRC)/IO/*.F90 $(BUILD)/.
	@cp $(SRC)/IO/*.c $(BUILD)/.
	@cp $(SRC)/Test_Suite/*.F90 $(BUILD)/.
	@cp $(SRC)/Physics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F $(BUILD)/.
	@cp $(SRC)/Include/*.F $(BUILD)/.
	@cp $(SRC)/Makefile $(BUILD)/.
	@cp $(SRC)/object_list $(BUILD)/.
ifeq ($(NODIRS),1)
	cp $(SRC)/Utility/MakeDir.F90_IBM $(BUILD)/MakeDir.F90
endif
ifdef CUSTOMROOT
	@echo Custom directory specified.
	@echo Any files from $(CUSTOMROOT) will overwrite standard Rayleigh source files.
	@cp $(CUSTOMROOT)/* $(BUILD)/. 2>/dev/null || :
endif

interp3d.gnu:
	@$(MAKE) --no-print-directory --directory=$(INTERP) interp3d.gnu
	@mkdir -p $(PREFIX)/bin
	@cp $(INTERP)/interp3d $(PREFIX)/bin/.

interp3d.intel:
	@$(MAKE) --no-print-directory --directory=$(INTERP) interp3d.intel
	@mkdir -p $(PREFIX)/bin
	@cp $(INTERP)/interp3d $(PREFIX)/bin/.


clean:
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean

clear_ipynb:
	etc/check_ipynb_cleared -c $^

.PHONY: install
install:
ifeq ($(target), "all")
	@echo "Installing executables into: " $(PREFIX)"/bin"
	@mkdir -p $(PREFIX)/bin
	@cp $(BUILD)/compiled/rayleigh.* $(PREFIX)/bin/.
else
ifdef output
	@echo "Installing rayleigh.$(target) into: " $(PREFIX)"/bin/$(output)"
	@mkdir -p $(PREFIX)/bin
	@cp $(BUILD)/compiled/rayleigh.$(target) $(PREFIX)/bin/$(output)
else
	@echo "Installing executables into: " $(PREFIX)"/bin"
	@mkdir -p $(PREFIX)/bin
	@cp $(BUILD)/compiled/rayleigh.* $(PREFIX)/bin/.
endif
endif

.PHONY: doc
doc:
	@sphinx-build -M html "." "doc/build"
	@sphinx-build -M latexpdf "." "doc/build"

.PHONY: distclean
distclean:
	rm -f $(BUILD)/compiled/rayleigh.*
	rm -f $(BUILD)/*.F
	rm -f $(BUILD)/*.F90
	rm -f $(BUILD)/*.c
	rm -f $(BUILD)/*.py
	rm -f $(BUILD)/*.o
	rm -f $(BUILD)/*.mod
	rm -f $(BUILD)/object_list
	rm -f $(BUILD)/Makefile
	rm -f $(BUILD)/Makefile.fdeps
	rm -f $(BUILD)/Machine_Definitions
	rm -f $(BUILD)/machine.no_comments
	@rm -f $(PREFIX)/bin/rayleigh.*
	@rm -f make.inc
	@rm -rf doc/build
	@echo "#Following configure, this file contains the definition of the PREFIX variable" >> make.inc


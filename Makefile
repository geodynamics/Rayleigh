include make.inc
SRC=src
BUILD=$(SRC)/build

# make the CUSTOMROOT variable available to sub-make processes
export CUSTOMROOT

rayleigh:
	@mkdir -p $(BUILD)/compiled
	@cp $(SRC)/Parallel_Framework/*.F90 $(BUILD)/.
	@cp $(SRC)/Data_Structures/*.F90 $(BUILD)/.
	@cp $(SRC)/Math_Layer/*.F90 $(BUILD)/.
	@cp $(SRC)/Plugins/*.F90 $(BUILD)/.
	@cp $(SRC)/IO/*.F90 $(BUILD)/.
	@cp $(SRC)/Test_Suite/*.F90 $(BUILD)/.
	@cp $(SRC)/Physics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F $(BUILD)/.
	@cp $(SRC)/Utility/*.F90 $(BUILD)/.
	@cp $(SRC)/Utility/*.c $(BUILD)/.
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
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) all
	@echo ""
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Compilation is complete."
	@echo " "
	@echo " "
	@echo "Run 'make install' to install the rayleigh executables into: "
	@echo $(PREFIX)"/bin"
	@echo ""
clean:
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean

.PHONY: install
install:
	@echo "Installing executables into: " $(PREFIX)"/bin"
	@mkdir -p $(PREFIX)/bin
	@cp $(BUILD)/compiled/rayleigh.* $(PREFIX)/bin/.

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
	rm -f $(BUILD)/Machine_Definitions
	rm -f $(BUILD)/machine.no_comments
	@rm -f $(PREFIX)/bin/rayleigh.*
	@rm -f make.inc
	@echo "#Following configure, this file contains the definition of the PREFIX variable" >> make.inc

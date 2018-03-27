include make.inc
SRC=src
BUILD=$(SRC)/build
rayleigh:
	@mkdir -p $(BUILD)/compiled
	@cp $(SRC)/parallel_framework/*.F90 $(BUILD)/.
	@cp $(SRC)/data_structures/*.F90 $(BUILD)/.
	@cp $(SRC)/math_layer/*.F90 $(BUILD)/.
	@cp $(SRC)/plugins/*.F90 $(BUILD)/.
	@cp $(SRC)/IO/*.F90 $(BUILD)/.
	@cp $(SRC)/test_suite/*.F90 $(BUILD)/.
	@cp $(SRC)/physics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F90 $(BUILD)/.
	@cp $(SRC)/Diagnostics/*.F $(BUILD)/.
	@cp $(SRC)/Utility/*.F90 $(BUILD)/.
	@cp $(SRC)/Utility/*.c $(BUILD)/.
	@cp $(SRC)/Makefile $(BUILD)/.
	@cp $(SRC)/object_list $(BUILD)/.
ifeq ($(NODIRS),1)
	cp $(SRC)/Utility/MakeDir.F90_IBM $(BUILD)/MakeDir.F90
endif
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) all
	@echo ""
	@echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	@echo "Compilation is complete."
	@echo "Run 'make install' to install the rayleigh executables into: "
	@echo $(PREFIX)"/bin"
	@echo ""
clean:
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean_exec
	@$(MAKE) --no-print-directory --directory=$(BUILD) clean


install:
	@echo "Installing executables into: " $(PREFIX)"/bin"
	@mkdir -p $(PREFIX)/bin
	@cp $(BUILD)/compiled/rayleigh.* $(PREFIX)/bin/.

distclean:
	rm -f $(BUILD)/compiled/rayleigh.*
	rm -f $(BUILD)/*.F
	rm -f $(BUILD)/*.F90
	rm -f $(BUILD)/*.c
	rm -f $(BUILD)/*.o
	rm -f $(BUILD)/*.mod
	rm -f $(BUILD)/Machine_Definitions
	rm -f $(BUILD)/machine.no_comments
	@rm -f $(PREFIX)/rayleigh.*
	@rm -f make.inc
	@echo "#Following configure, this file contains the definition of the PREFIX variable" >> make.inc

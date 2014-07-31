##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

CMAKE ?= cmake

PFUNIT_SOURCE = $(abspath ../pfunit)
PFUNIT_BUILD = $(BUILD_DIR)/pfunit
export PFUNIT_INSTALL = $(BUILD_DIR)/pfunit-install

ifeq '$(FC)' 'ifort'
  PFUNIT_COMPILER_ID = Intel
else ifeq '$(FC)' 'gfortran'
  PFUNIT_COMPILER_ID = GNU
else ifeq '$(FC)' 'nagfor'
  PFUNIT_COMPILER_ID = NAG
else ifeq '$(FC)' 'xlf'
  PFUNIT_COMPILER_ID = XL
endif

.PHONY: pfunit
pfunit: compilertest $(PFUNIT_BUILD)
	$(MAKE) -C $(PFUNIT_BUILD) install

$(PFUNIT_BUILD):
	mkdir $@
	cd $@; $(CMAKE) -DCMAKE_Fortran_COMPILER=$(FC) \
	                -DCMAKE_INSTALL_PREFIX=$(PFUNIT_INSTALL) \
	                $(PFUNIT_SOURCE)

.PHONY: compilertest
ifdef IFORT_VERSION
compilertest:
	@if [ $(IFORT_VERSION) -lt 0130000 ]; then \
	  echo >&2 "*** [ERROR] pFUnit will only compile with ifort v13 or later."; \
	  false; \
	fi
else ifdef GFORTRAN_VERSION
compilertest:
	echo $(GFORTRAN_VERSION)
	@if [ $(GFORTRAN_VERSION) -lt 040500 ]; then \
	  echo >&2 "*** [ERROR] pFUnit will only compile with gfortran v4.5 or later."; \
	  false; \
	fi
else
compilertest:
endif

.PHONY: cleanpfunit
cleanpfunit:
	-rm -rf $(PFUNIT_BUILD)
	-rm -rf $(PFUNIT_INSTALL)

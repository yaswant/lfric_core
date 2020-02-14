##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = lfric_atm

PROFILE ?= fast-debug
LFRIC_TARGET_PLATFORM ?= generic-x86

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5 \
                                    $(CXX_RUNTIME_LIBRARY)
export EXTERNAL_STATIC_LIBRARIES = xios

export UMPHYSICS_TARGETS = extract-um-physics
export UMPHYSICS_CLEAN   = clean-um-physics
export UM_FCM=$(UM_PHYS_DIR)/
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..
export PATH        := $(PATH):$(abspath $(ROOT_DIR)/GPL/bin)
export PYTHONPATH  := $(PYTHONPATH):$(abspath $(ROOT_DIR)/GPL/lib/python)
export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include build/fortran/$(FORTRAN_COMPILER).mk
include build/project.mk

export FFLAGS  := $(FFLAGS) $(OPENMP_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)

##############################################################################
# Launch gscan to monitor suites from this suite-group
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan  $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &           \
          usleep 1                                                                ;\
        fi

##############################################################################
# Test
#
.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Build
#
.PHONY: build
ifeq "$(PROFILE)" "full-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_INIT) $(FFLAGS_RUNTIME) $(FFLAGS_NO_OPTIMISATION)
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_SAFE_OPTIMISATION)
else ifeq "$(PROFILE)" "production"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_RISKY_OPTIMISATION)
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
build: PROJECT      = lfric_atm
build: SCRATCH_DIR := $(WORKING_DIR)/um-scratch
build: WORKING_DIR := $(WORKING_DIR)/lfric_atm
build: export platform_config_dir=$(strip $(shell grep -E "$(LFRIC_TARGET_PLATFORM)\s*:\s*$(FORTRAN_COMPILER)\s*:" $(PROJECT_DIR)/fcm-make/target-map.txt | cut -d : -f 3))
build: export optimisation_level=$(strip $(shell grep -E "$(PROFILE)\s*:\s*" $(PROJECT_DIR)/fcm-make/optimisation-map.txt | cut -d : -f 2))
build: ALWAYS
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	         SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	         WORKING_DIR=$(WORKING_DIR) \
	         WITHOUT_PROGRAMS=1
	$(call MESSAGE,====================================)
	$(call MESSAGE,Psycloning Gungho dynamical core)
	$(call MESSAGE,====================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	        SOURCE_DIR=$(ROOT_DIR)/gungho/source \
                WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=gungho \
                  SOURCE_DIR=$(ROOT_DIR)/gungho/source \
                  WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$(call MESSAGE,Extracting,"UM physics" )
	$Qif [ x$(platform_config_dir) = x ] ; then echo Unable to convert $(LFRIC_TARGET_PLATFORM):$(FORTRAN_COMPILER) to a UM target; false; fi
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source \
	          WORKING_DIR=$(WORKING_DIR)
	$(info Using UM target $(platform_config_dir))
	# Retrieve and preprocess the UM, Jules and Socrates code
	# The UM_ENV file contains the appropriate locations and UM side
	# environment variables
	$Qsource $(PROJECT_DIR)/fcm-make/parameters.sh && \
	    fcm make -C $(SCRATCH_DIR) -f $(PROJECT_DIR)/fcm-make/extract.cfg
	# Note that if wanting to modify UM source this should be done via the
	# UM repository either through a working copy or branch
	$Qrsync -acvz $(SCRATCH_DIR)/preprocess-atmos/ $(WORKING_DIR)/science/
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
                  BIN_DIR=$(BIN_DIR) \
                  PROGRAMS="$(PROGRAMS)"


##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export ADDITIONAL_EXTRACTION = $(ROOT_DIR)/infrastructure/source
unit-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests: export PROGRAMS = lfric_atm_unit_tests
unit-tests: export PROJECT = lfric_atm 
unit-tests: export SOURCE_DIR = source
unit-tests: export TEST_DIR = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: do-unit-tests
	$(Q)echo > /dev/null


##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export ADDITIONAL_EXTRACTION = $(ROOT_DIR)/infrastructure/source
integration-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
integration-tests: export PROJECT =lfric_atm 
integration-tests: export SOURCE_DIR = source
integration-tests: export TEST_DIR = unit-test
integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-tests
integration-tests: do-integration-tests
	$(Q)echo > /dev/null


##############################################################################
# Clean
#
.PHONY: clean
clean: partial-clean
	$(call MESSAGE,Removing,"UM work space")
	$(Q)-rm -r $(WORKING_DIR)

.PHONY: partial-clean pclean
partial-clean pclean: ALWAYS
	$(call MESSAGE,Removing,"LFRic atmosphere  workspace")
	$(Q)-rm -r $(WORKING_DIR)/lfric_atm
	$(call MESSAGE,Removing,"LFRic atmosphere binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi


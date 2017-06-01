##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Root make file for GungHo. Targets provided our detailed below...
#
# all: (default) Complete build and test of GungHo.
#      Intended for interactive use.
# build: Build all executables and any supporting files
# test: Run test battery including unit tests and others
# clean: Delete all final products and working files
#
# The following variables may be specified to modify the build process...
#
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PROFILE: Set to a string representing a package of compiler options.
#          Potential profiles are 'full-debug', 'fast-debug' and 'production'.
#          Default: fast-debug
# LINK_TYPE: Either 'static' or 'dynamic'.
#            Default: dynamic
# OPTIMISATION_PROFILE: PSyclone will use optimisation scripts for this
#                       platform.
#
##############################################################################

export IGNORE_DEPENDENCIES = netcdf MPI ESMF pfunit_mod
export EXTERNAL_DYNAMIC_LIBRARIES = esmf netcdff netcdf hdf5

export ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

ifdef OPTIMISATION_PROFILE
  export OPTIMISATION_PATH = gungho/optimisation/$(OPTIMISATION_PROFILE)
endif

PROFILE ?= fast-debug

.PHONY: all
all:
	$(MAKE) build
	$(MAKE) test

.PHONY: build
build: build-gungho build-mesh_tools

.PHONY: test
test: test-infrastructure test-gungho test-mesh_tools

.PHONY: test-suite
test-suite: SUITE_GROUP ?= developer
test-suite:
	$(Q)if [ -z "$(TEST_SUITE_TARGETS)" ] ; then \
	    echo *** Please set the DYNAMO_TEST_SUITE_TARGETS environment variable. ; \
	    exit 1 ; \
	fi
	$(Q)umask 022; for target in $(DYNAMO_TEST_SUITE_TARGETS) ; do \
	    echo Launching test suite against $$target ; \
	    rose stem --name=$(shell basename `pwd`)-infrastructure-$$target-$(SUITE_GROUP) --config=infrastructure/rose-stem --opt-conf-key=$$target --group=$(SUITE_GROUP); \
	    rose stem --name=$(shell basename `pwd`)-gungho-$$target-$(SUITE_GROUP) --config=gungho/rose-stem --opt-conf-key=$$target --group=$(SUITE_GROUP); \
	    rose stem --name=$(shell basename `pwd`)-mesh_tools-$$target-$(SUITE_GROUP) --config=mesh_tools/rose-stem --opt-conf-key=$$target --group=$(SUITE_GROUP); \
	done

include $(ROOT)/infrastructure/build/lfric.mk

export FFLAGS  := $(FFLAGS) $(OPENMP_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)

##############################################################################
# Infrastructure
#
.PHONY: extract-infrastructure
extract-infrastructure:export SOURCE_DIR   = infrastructure/source
extract-infrastructure:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Test infrastructure
#
.PHONY: test-infrastructure
test-infrastructure: unit-test-infrastructure integration-test-infrastructure

.PHONY: integration-test-infrastructure
integration-test-infrastructure: export PROJECT      = infrastructure
integration-test-infrastructure: export SOURCE_DIR   = infrastructure/integration-test
integration-test-infrastructure: export WORKING_DIR := $(WORKING_DIR)/infrastructure
integration-test-infrastructure: export BIN_DIR      = $(ROOT)/tests
integration-test-infrastructure: export PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
integration-test-infrastructure: generate-integration-test-infrastructure
	$(MAKE) run-integration-test-infrastructure

.PHONY: generate-integration-test-infrastructure
generate-integration-test-infrastructure: extract-integration-test-infrastructure \
                                         extract-infrastructure
	$(Q)echo >/dev/null

.PHONY: extract-integration-test-infrastructure
extract-integration-test-infrastructure:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk

.PHONY: unit-tests-infrastructure
unit-test-infrastructure: export EXTERNAL_STATIC_LIBRARIES = pfunit
unit-test-infrastructure: export PROJECT      = infrastructure
unit-test-infrastructure: export WORKING_DIR := $(WORKING_DIR)/infrastructure
unit-test-infrastructure: export BIN_DIR      = $(ROOT)/tests
unit-test-infrastructure: export PROGRAMS     = infrastructure_unit_test
unit-test-infrastructure: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-test-infrastructure: export PRE_PROCESS_MACROS = \
                         PFUNIT_EXTRA_USAGE=infrastructure_suite_fixture_mod \
                         PFUNIT_EXTRA_INITIALIZE=set_up_suite \
                         PFUNIT_EXTRA_FINALIZE=tear_down_suite
unit-test-infrastructure: generate-unit-test-infrastructure
	$(MAKE) run-unit-test-infrastructure

.PHONY: generate-unit-test-infrastructure
generate-unit-test-infrastructure: export SOURCE_DIR = infrastructure/unit-test
generate-unit-test-infrastructure: extract-unit-test-infrastructure \
                                   extract-infrastructure
	$(MAKE) pfunit

.PHONY: extract-unit-test-infrastructure
extract-unit-test-infrastructure: export SOURCE_DIR = infrastructure/unit-test
extract-unit-test-infrastructure:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk
	$(Q) cp -r $(SOURCE_DIR)/data $(WORKING_DIR)/

##############################################################################
# Build GungHo
#
.PHONY: build-gungho
ifeq "$(PROFILE)" "full-debug"
build-gungho: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_INIT) $(FFLAGS_RUNTIME) $(FFLAGS_NO_OPTIMISATION)
else ifeq "$(PROFILE)" "fast-debug"
build-gungho: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_SAFE_OPTIMISATION)
else ifeq "$(PROFILE)" "production"
build-gungho: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_RISKY_OPTIMISATION)
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build-gungho: export PROJECT      = gungho
build-gungho: export SOURCE_DIR   = gungho/source
build-gungho: export WORKING_DIR := $(WORKING_DIR)/gungho
build-gungho: export DATABASE     = $(abspath $(WORKING_DIR)/dependencies.db)
build-gungho: export BIN_DIR      = $(ROOT)/bin
build-gungho: export PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
build-gungho: generate-configuration-gungho generate-psykal-gungho
	$(MAKE) -f $(LFRIC_BUILD)/lfric.mk compile

.PHONY: generate-configuration-gungho
generate-configuration-gungho: export SOURCE_DIR = gungho/source
generate-configuration-gungho: extract-gungho extract-infrastructure
	$(MAKE) -f $(LFRIC_BUILD)/lfric.mk configuration

.PHONY: generate-psykal-gungho
generate-psykal-gungho: export SOURCE_DIR   = gungho/source
generate-psykal-gungho: extract-gungho extract-infrastructure
	$(MAKE) -f $(LFRIC_BUILD)/lfric.mk psykal

.PHONY: extract-gungho
extract-gungho: export SOURCE_DIR = gungho/source
extract-gungho:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Test GungHo
#
.PHONY: test-gungho
test-gungho: unit-test-gungho

.PHONY: unit-tests-gungho
unit-test-gungho: export PROJECT      = gungho
unit-test-gungho: export WORKING_DIR := $(WORKING_DIR)/gungho
unit-test-gungho: export BIN_DIR      = $(ROOT)/tests
unit-test-gungho: export PROGRAMS     = gungho_unit_test
unit-test-gungho: export EXTERNAL_STATIC_LIBRARIES = pfunit
unit-test-gungho: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-test-gungho: export PRE_PROCESS_MACROS = \
                                 PFUNIT_EXTRA_USAGE=gungho_suite_fixture_mod \
                                 PFUNIT_EXTRA_INITIALIZE=set_up_suite \
                                 PFUNIT_EXTRA_FINALIZE=tear_down_suite
unit-test-gungho: generate-configuration-gungho generate-psykal-gungho \
                  generate-unit-test-gungho
	$(MAKE) run-unit-test-gungho

.PHONY: generate-unit-test-gungho
generate-unit-test-gungho: export SOURCE_DIR   = gungho/unit-test
generate-unit-test-gungho: extract-unit-test-gungho \
                           extract-gungho extract-infrastructure
	$(MAKE) pfunit

.PHONY: extract-unit-test-gungho
extract-unit-test-gungho: export SOURCE_DIR = gungho/unit-test
extract-unit-test-gungho:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk
	$(Q)cp -r $(SOURCE_DIR)/data $(WORKING_DIR)

##############################################################################
# Build mesh_tools
#
.PHONY: build-mesh_tools
ifeq "$(PROFILE)" "full-debug"
build-mesh_tools: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_INIT) $(FFLAGS_RUNTIME) $(FFLAGS_NO_OPTIMISATION)
else ifeq "$(PROFILE)" "fast-debug"
build-mesh_tools: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_SAFE_OPTIMISATION)
else ifeq "$(PROFILE)" "production"
build-mesh_tools: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_RISKY_OPTIMISATION)
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build-mesh_tools: export PROJECT      = mesh_tools
build-mesh_tools: export SOURCE_DIR   = mesh_tools/source
build-mesh_tools: export WORKING_DIR := $(WORKING_DIR)/mesh_tools
build-mesh_tools: export DATABASE     = $(abspath $(WORKING_DIR)/dependencies.db)
build-mesh_tools: export BIN_DIR      = $(ROOT)/bin
build-mesh_tools: export PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
build-mesh_tools: generate-configuration-mesh_tools
	$(MAKE) -f $(LFRIC_BUILD)/lfric.mk compile

.PHONY: generate-configuration-mesh_tools
generate-configuration-mesh_tools: export SOURCE_DIR   = mesh_tools/source
generate-configuration-mesh_tools: extract-mesh_tools
	$(MAKE) extract-infrastructure
	$(MAKE) -f $(LFRIC_BUILD)/lfric.mk configuration

.PHONY: extract-mesh_tools
extract-mesh_tools: export SOURCE_DIR   = mesh_tools/source
extract-mesh_tools: export WORKING_DIR := $(WORKING_DIR)/mesh_tools
extract-mesh_tools:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Test mesh_tools
#
.PHONY: test-mesh_tools
test-mesh_tools: export PROJECT=mesh_tools
test-mesh_tools: unit-test-mesh_tools

.PHONY: unit-tests-mesh_tools
unit-test-mesh_tools: export EXTERNAL_STATIC_LIBRARIES = pfunit
unit-test-mesh_tools: export WORKING_DIR := $(WORKING_DIR)/mesh_tools
unit-test-mesh_tools: export BIN_DIR      = $(ROOT)/tests
unit-test-mesh_tools: export PROGRAMS     = mesh_tools_unit_test
unit-test-mesh_tools: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-test-mesh_tools: export PRE_PROCESS_MACROS = \
                                PFUNIT_EXTRA_USAGE=mesh_tools_suite_fixture_mod \
                                PFUNIT_EXTRA_INITIALIZE=set_up_suite \
                                PFUNIT_EXTRA_FINALIZE=tear_down_suite
unit-test-mesh_tools: generate-configuration-mesh_tools generate-unit-test-mesh_tools
	$(MAKE) run-unit-test-mesh_tools

.PHONY: generate-unit-test-mesh_tools
generate-unit-test-mesh_tools: export SOURCE_DIR   = mesh_tools/unit-test
generate-unit-test-mesh_tools: extract-unit-test-mesh_tools \
                            extract-mesh_tools extract-infrastructure
	$(MAKE) pfunit

.PHONY: extract-unit-test-mesh_tools
extract-unit-test-mesh_tools: export SOURCE_DIR = mesh_tools/unit-test
extract-unit-test-mesh_tools:
	$(MAKE) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Remove working files
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,work space)
	$(Q)-rm -rf $(WORKING_DIR)
	$(call MESSAGE,Removing,binaries)
	$(Q)-rm -rf bin
	$(call MESSAGE,Removing,tests)
	$(Q)-rm -rf tests

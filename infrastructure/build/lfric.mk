##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Include this file from your model make file in order to gain access to the
# LFRic build system. Include it at the end of the make file as it contains
# targets which you do not want to become the default target.
#
# Variables provided by including this file...
#
# LFRIC_BUILD: Path to the build system
# COMPILE_OPTIONS: File of target-specific compile options used in compile.mk
#
# Macros expected by the build system are as follows...
#
# ROOT_DIR: Absolute path to the project directory
# OPTIMISATION_PATH: Where PSyclone optimisation scripts may be found.
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PURGE_SUITES: Set to non-zero value to clean out exisiting rose suites of
#               same name. (this is also the default action)
#               Set to zero to not clean out existing rose suite.
# TEST_SUITE_TARGETS: Space separated list of target identifiers to be used
#                     when launching the test suite. Default is "meto-spice
#                     meto-xc40"
#
##############################################################################

.SECONDEXPANSION:

# Ensure make offers the features we need...
#
$(info ** Make version $(MAKE_VERSION))
ifeq ($(filter else-if,$(value .FEATURES)),)
  $(error The build system requires else-if support from GMake)
endif

# Default variables...
#
export WORKING_DIR ?= working
export PWD ?= $(shell pwd)

TEST_SUITE_TARGETS ?= meto-spice meto-xc40

# Make the build system available...
#
export LFRIC_BUILD := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# Make the infrastructure available...
#
export LFRIC_INFRASTRUCTURE := $(realpath $(LFRIC_BUILD)/..)

# Attempt to identify Cray systems...
#
ifdef PE_ENV
  CRAY_ENVIRONMENT = true
  export CRAY_ENVIRONMENT
endif

# Set the default precision for reals
RDEF_PRECISION ?= 64
export PRE_PROCESS_MACROS += RDEF_PRECISION=$(RDEF_PRECISION)

# The compile options file overrides compile options based on file-name pattern matching.
# Use the miniapp-specific file if it exists. Otherwise use the infrastructure file.
ifeq (,$(wildcard $(PROJECT_DIR)/build/compile_options.mk))
  export COMPILE_OPTIONS := $(abspath $(ROOT_DIR)/infrastructure/build/compile_options.mk)
else
  export COMPILE_OPTIONS := $(abspath $(PROJECT_DIR)/build/compile_options.mk)
endif

# Set up verbose logging...
#
ifdef VERBOSE
  Q :=
  VERBOSE_ARG = -verbose
  SHORT_VERBOSE_ARG = -v
  DOUBLE_VERBOSE_ARG = --verbose
else
  Q := @
  QUIET_ARG = --quiet
  VERBOSE_REDIRECT = >/dev/null
endif

export Q QUIET_ARG VERBOSE_REDIRECT

# Set flag to perform a fresh rose stem suite

CLEAN_OPT ?= '--new'
ifeq '$(PURGE_SUITES)' '0'
  CLEAN_OPT =
endif

# Fix enumeration arguments between builds, only for technical
# development

FIX_ENUMS ?= 0
ifeq '$(FIX_ENUMS)' '1'
  FIX_ENUMS_OPT = '--norandom_enums'
else
  FIX_ENUMS_OPT =
endif

# We only want to send terminal control characters if there is a terminal to
# interpret them...
#
_NOW = `date +%H:%M:%S`
ifneq 'x$(TERM)' 'x'
  MESSAGE = $(Q)printf "%s \\033[1m$(1)\\033[0m %s\n" $(_NOW) $(2)
else
  MESSAGE = $(Q)echo $(_NOW) *$(1)* $(2)
endif

# Set up some special macros for hard to escape characters
#
EMPTY :=
SPACE := $(EMPTY) # This comment highlights space character.
PERCENT := %
OPEN_PAREN := (

# Prerequisite for targets which should always be run.
#
.PHONY: ALWAYS
ALWAYS:

# The directory containing the target. Useful for order-only prerequisites to
# create that directory.
#
TARGET_DIR = $(patsubst $(PERCENT)/,$(PERCENT),$(dir $@))

##############################################################################
# Build UML documentation
#
# DOCUMENT_DIR - Directory in which documentation is placed.
# SOURCE_DIR   - Root directory of source tree.
# WORKING_DIR  - Location for temporary files.
#
.PHONY: uml-documentation
uml-documentation:
	$(Q)$(MAKE) $(QUIET_ARG) uml-pdfs DOCUMENT_DIR=$(DOCUMENT_DIR) SOURCE_DIR=$(SOURCE_DIR) WORKING_DIR=$(WORKING_DIR)

.PHONY: uml-pdfs
uml-pdfs: $$(patsubst $$(SOURCE_DIR)/$$(PERCENT).puml, \
                      $$(DOCUMENT_DIR)/$$(PERCENT).pdf, \
                      $$(wildcard $$(SOURCE_DIR)/*.puml))
	$(Q)echo >/dev/null

.PRECIOUS: $(DOCUMENT_DIR)/%.pdf
$(DOCUMENT_DIR)/%.pdf: $(DOCUMENT_DIR)/%.svg
	$(call MESSAGE,Translating,$(notdir $<))
	$Qconvert $< $@

.PRECIOUS: $(DOCUMENT_DIR)/%.svg
$(DOCUMENT_DIR)/%.svg: $(SOURCE_DIR)/%.puml \
                       $$(addprefix $$(SOURCE_DIR)/, \
                                    $$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$(SOURCE_DIR)/$$*.puml))
	$(call MESSAGE,Generating,$(notdir $@))
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)plantuml $(SHORT_VERBOSE_ARG) -tsvg -o $(abspath $(dir $@)) $(abspath $<)

##############################################################################
# Build API documentation
#
# PROJECT      - Name of the project for logging purposes.
# DOCUMENT_DIR - Directory in which documentation is placed.
# CONFIG_DIR   - Directory where Doxyfile can be found.
# SOURCE_DIR   - Root directory of source tree.
# WORKING_DIR  - Location for temporary files.
#
.PHONY: api-documentation
api-documentation: ALWAYS
	$(call MESSAGE,API,$(PROJECT))
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)( cat $(CONFIG_DIR)/Doxyfile; \
	      echo INPUT=$(SOURCE_DIR) $(CONFIG_DIR)/$$(sed -n -e 's/\s*INPUT\s*=\s*//p' $(CONFIG_DIR)/Doxyfile); \
	      echo USE_MDFILE_AS_MAINPAGE=$(CONFIG_DIR)/$$(sed -n -e 's/\s*USE_MDFILE_AS_MAINPAGE\s*=\s*//p' $(CONFIG_DIR)/Doxyfile); \
	      echo OUTPUT_DIRECTORY=$(DOCUMENT_DIR) ) \
	    | doxygen - $(VERBOSE_REDIRECT)

##############################################################################
# Launch test suite
#
# SUITE_CONFIG    - Path to rose-stem directory.
# SUITE_BASE_NAME - Name for suites.
#
.PHONY: launch-test-suite
launch-test-suite: SUITE_GROUP ?= developer
launch-test-suite: SUITE_NAME = $(SUITE_BASE_NAME)-$$target-$(SUITE_GROUP)
ifdef VERBOSE
launch-test-suite: VERBOSE_ARG = --define-suite=VERBOSE=$(VERBOSE)
endif
launch-test-suite:
	$(Q)umask 022; for target in $(TEST_SUITE_TARGETS) ; do \
	echo Launching $(PROJECT_NAME) test suite against $$target with $(SUITE_GROUP) group ; \
	rose stem --name=$(SUITE_NAME) \
	          --config=$(SUITE_CONFIG) \
	          --opt-conf-key=$$target --no-gcontrol \
	          $(CLEAN_OPT) $(QUIET_ARG) \
	          --define-suite=RDEF_PRECISION=$(RDEF_PRECISION) \
                  $(VERBOSE_ARG) \
	          --group=$(SUITE_GROUP); \
	done

##############################################################################
# Extract parts of other projects.
#
.PHONY: %/extract
%/extract:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	            SOURCE_DIR=$* WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Invoke PSyclone to generate PSy layer.
#
# Psyclone is called on the original source but that source may use other
# modules so extraction must be complete first.
#
.PHONY: %/psyclone
%/psyclone: $$(addsuffix /extract, $$*)
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	            SOURCE_DIR=$* \
	            WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Invoke Configurator to generate namelist handling code.
#
# Configurator is called on the original source but that source may use other
# modules so extraction must be complete first.
#
.PHONY: %/configure
%/configure: PROJECT ?= $(lastword $(filter-out source, $(subst /, ,$*)))
%/configure: $$(addsuffix /extract, $$*)
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	            PROJECT=$(PROJECT) \
	            SOURCE_DIR=$* \
	            WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Run unit tests.
#
.PHONY: do-unit-tests
do-unit-tests: export ANY_TESTS := $(strip $(shell find $(SOURCE_DIR) -name *test.pf))
do-unit-tests:
	$(Q)$(MAKE) $(QUIET_ARG) do-unit-test/detect \
	            ADDITIONAL_EXTRACTION=$(ADDITIONAL_EXTRACTION) \
	            ANY_TESTS=$(ANY_TESTS) \
	            BIN_DIR=$(BIN_DIR) \
	            EXTERNAL_STATIC_LIBRARIES="$(EXTERNAL_STATIC_LIBRARIES)" \
	            PROGRAMS=$(PROGRAMS) \
	            PROJECT=$(PROJECT) \
	            SOURCE_DIR=$(SOURCE_DIR) \
	            TEST_DIR=$(TEST_DIR) \
	            WORKING_DIR=$(WORKING_DIR)

.PHONY: do-unit-test/%
do-unit-test/detect: $(if $(ANY_TESTS), do-unit-test/found, \
                                        do-unit-test/none)
	$(Q)echo >/dev/null

do-unit-test/none:
	$(Q)$(call MESSAGE,Unit tests,'None to run')

do-unit-test/found: do-unit-test/build
	$(call MESSAGE,Running,$(PROGRAMS))
	$(Q)cd $(WORKING_DIR); \
            mpiexec -n 4 $(BIN_DIR)/$(PROGRAMS) $(DOUBLE_VERBOSE_ARG)

do-unit-test/build: WITHOUT_PROGRAMS = 1
do-unit-test/build: PROCESS_TARGETS = $(SOURCE_DIR) $(ADDITIONAL_EXTRACTION)
do-unit-test/build: $$(addsuffix /configure, $$(PROCESS_TARGETS)) \
                    $$(addsuffix /psyclone, $$(PROCESS_TARGETS)) \
                    $$(addsuffix /extract, $$(TEST_DIR))
	$(Q)mkdir -p $(WORKING_DIR)
	$(Q)if [ -d $(TEST_DIR)/data ] ; then \
	        cp -r $(TEST_DIR)/data $(WORKING_DIR); fi
	$(Q)$(MAKE) $(QUIET_ARG) pfunit \
	            SOURCE_DIR=$(TEST_DIR) WORKING_DIR=$(WORKING_DIR)
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(Q)$(MAKE) $(QUIET_ARG) \
                    -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
                    PRE_PROCESS_MACROS=$(PRE_PROCESS_MACROS)


##############################################################################
# Run integration tests.
#
.PHONY: do-integration-tests
do-integration-tests: export ANY_TESTS = $(patsubst $(TEST_DIR)/%,%,$(basename $(shell find $(TEST_DIR) -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \; 2>/dev/null)))
do-integration-tests:
	$(Q)$(MAKE) $(QUIET_ARG) do-integration-test/detect \
	            ADDITIONAL_EXTRACTION=$(ADDITIONAL_EXTRACTION) \
	            ANY_TESTS="$(ANY_TESTS)" \
	            BIN_DIR=$(BIN_DIR) \
	            META_FILE_DIR=$(META_FILE_DIR) \
	            PRE_PROCESS_MACROS=$(PRE_PROCESS_MACROS) \
	            PROGRAMS="$(ANY_TESTS)" \
	            PROJECT_NAME=$(PROJECT_NAME) \
	            SOURCE_DIR=$(SOURCE_DIR) \
	            TEST_DIR=$(TEST_DIR) \
	            WORKING_DIR=$(WORKING_DIR)

.PHONY: do-integration-test/%
do-integration-test/detect: $(if $(ANY_TESTS), \
                                 $(ANY_TESTS:%=do-integration-test/run/%), \
                                 do-integration-test/none)
	$(Q)echo >/dev/null

do-integration-test/none:
	$(Q)$(call MESSAGE,Integration tests,'None to run')

do-integration-test/run/%: export PYTHONPATH := $(PYTHONPATH):$(LFRIC_BUILD)
do-integration-test/run/%: do-integration-test/build
	$(call MESSAGE,Running,$*)
	$(Q)cd $(TEST_DIR)/$(dir $*); \
	    ./$(notdir $(addsuffix .py,$*)) $(addprefix $(BIN_DIR)/,$(notdir $*))

do-integration-test/build: PROCESS_TARGETS = $(SOURCE_DIR) $(TEST_DIR) $(ADDITIONAL_EXTRACTION)
do-integration-test/build: $$(addsuffix /configure, $$(PROCESS_TARGETS)) \
                           $$(addsuffix /psyclone, $$(PROCESS_TARGETS))
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk


##############################################################################
# Generate configuration source.
#
.PHONY: configuration
configuration: FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
configuration:
	$(call MESSAGE,lfric.mk $(FIX_ENUMS_OPT))
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk

##############################################################################
# Generate pFUnit unit tests.
#
.PHONY: pfunit
pfunit:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk

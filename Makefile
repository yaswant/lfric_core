##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

# By default we link dynamically
#
LINK ?= DYNAMIC

# A trio of targets are offered for various levels of debugging support and
# speed of execution.
#
.PHONY: fast-debug
fast-debug: OPTIMISATION?=SAFE
fast-debug: SYMBOLS?=YES
fast-debug: test

.PHONY: full-debug
full-debug: OPTIMISATION?=NONE
full-debug: SYMBOLS?=YES
full-debug: CHECKS?=YES
full-debug: test

.PHONY: production
production: OPTIMISATION?=RISKY
production: SYMBOLS?=YES
production: test

# The 'build' target has some default behaviour with regard to build
# type. This is to support the situation where 'test' is called on
# its own. It also allows 'build' to be called on its own.
#
.PHONY: build
build: OPTIMISATION?=SAFE
build: SYMBOLS?=YES
build: tools scripts
	$(MAKE) -C src/dynamo LINK=$(LINK) OPTIMISATION=$(OPTIMISATION) \
	                      SYMBOLS=$(SYMBOLS) CHECKS=$(CHECKS)

.PHONY: tools
tools:
	$(MAKE) -C tools

.PHONY: scripts
scripts:
	$(MAKE) -C src/scripts

# The 'test' target allows tests to be run on their own. When it is used in
# this sense there is no way to know which (if any) build target was used.
# In this case the build default will be used.
#
.PHONY: test
test: build
ifndef NOTEST
	$(MAKE) -C src/test
	$(MAKE) -C src/functional-test
endif

.PHONY:test-suite
test-suite:
	@if [ -z "$(DYNAMO_TEST_SUITE_TARGETS)" ] ; then \
	    echo *** Please set the DYNAMO_TEST_SUITE_TARGETS environment variable. ; \
	    exit 1 ; \
	fi
	@for target in $(DYNAMO_TEST_SUITE_TARGETS) ; do \
	    echo Launching test suite against $$target ; \
	    rose stem --name=$(shell basename `pwd`)-$$target-d --opt-conf-key=$$target ; \
	done

.PHONY:nightly-suite
nightly-suite:
	@if [ -z "$(DYNAMO_TEST_SUITE_TARGETS)" ] ; then \
	    echo *** Please set the DYNAMO_TEST_SUITE_TARGETS environment variable. ; \
	    exit 1 ; \
	fi
	@for target in $(DYNAMO_TEST_SUITE_TARGETS) ; do \
	    echo Launching nightly test suite against $$target ; \
	    rose stem --name=$(shell basename `pwd`)-$$target-n --group=nightly --opt-conf-key=$$target ; \
	done

.PHONY:csar-suite
csar-suite:
	@if [ -z "$(DYNAMO_TEST_SUITE_TARGETS)" ] ; then \
	    echo *** Please set the DYNAMO_TEST_SUITE_TARGETS environment variable. ; \
	    exit 1 ; \
	fi
	@for target in $(DYNAMO_TEST_SUITE_TARGETS) ; do \
	    echo Launching csar test suite against $$target ; \
	    rose stem --name=$(shell basename `pwd`)-$$target-csar --group=csar --opt-conf-key=$$target ; \
	done

# Build the projects documentation. This includes both API and design documents.
.PHONY: doc docs
doc docs:
	$(MAKE) -C documentation

# Clean only the dynamo build. This leaves pFUnit alone.
#
.PHONY: clean
clean:
	$(MAKE) -C src/dynamo clean
	$(MAKE) -C src/test clean
	$(MAKE) -C src/functional-test clean

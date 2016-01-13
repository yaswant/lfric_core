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
build: tools
	$(MAKE) -C src/dynamo LINK=$(LINK) OPTIMISATION=$(OPTIMISATION) \
	                      SYMBOLS=$(SYMBOLS) CHECKS=$(CHECKS)

.PHONY: tools
tools:
	$(MAKE) -C tools

# The 'test' target allows tests to be run on their own. When it is used in
# this sense there is no way to know which (if any) build target was used.
# In this case the build default will be used.
#
.PHONY: test
test: build
	$(MAKE) -C src/test

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

# Clean both dynamo and pFUnit builds.
#
.PHONY: clean-all
clean-all:
	$(MAKE) -C src/dynamo clean ALL=1
	$(MAKE) -C src/test clean ALL=1

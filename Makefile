##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# LFRic project make file. This simply calls down into sub-project make files.
#
# Variables:
#   OPERATE_ON - Sub-projects which will be affected by operations.
#                default: infrastructure, mesh_tools and gungho
#   TEST_SUITE_TARGETS - Platforms to target with test suite.
#
#############################################################################

# Operate only on this list of sub-projects. May be overridden from the
# terminal.
#
OPERATE_ON ?= infrastructure mesh_tools gungho lfric_atm \
              miniapps/skeleton                          \
              miniapps/gravity_wave

##############################################################################
# Perform default action on each sub-project in OPERATE_ON list.
#
.PHONY: default
default: $(addprefix default/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: default/%
default/%: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -C $*

##############################################################################
# Perform the clean action on each sub-project in OPERATE_ON list.
#
.PHONY: clean
clean: $(addprefix clean/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: clean/%
clean/%: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -C $* clean

##############################################################################
# Launch test suite for each sub-project in OPERATE_ON list.
#
.PHONY: test-suite
test-suite: SUITE_GROUP ?= developer
test-suite: $(addprefix test-suite/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: test-suite/%
test-suite/%: ALWAYS
	$(Q)-$(MAKE) $(QUIET_ARG) -C $* test-suite TEST_SUITE_TARGETS="$(TEST_SUITE_TARGETS)"

##############################################################################

include infrastructure/build/lfric.mk

##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

MPIC_COMPILER := $(shell $(CXX) --version \
                   | awk -F " " 'NR<2 { printf "%s", $$1 }')
$(info ** Chosen MPI C++ compiler "$(MPIC_COMPILER)")

ifeq '$(MPIC_COMPILER)' 'g++'
  CXX_COMPILER = g++
else ifeq '$(MPIC_COMPILER)' 'icc'
  CXX_COMPILER = icc
else ifeq '$(MPIC_COMPILER)' 'Cray'
  CXX_COMPILER = craycc
else
  $(error Unrecognised mpic++ compiler option: "$MPIC_COMPILER")
endif

include $(LFRIC_BUILD)/cxx/$(CXX_COMPILER).mk

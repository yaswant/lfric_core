##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

MPIF90_COMPILER := $(shell $(FC) --version \
                   | awk -F " " 'NR<2 { printf "%s", $$1 }')
$(info ** Chosen MPI Fortran compiler "$(MPIF90_COMPILER)")

ifeq '$(MPIF90_COMPILER)' 'GNU'
  FORTRAN_COMPILER = gfortran
else ifeq '$(MPIF90_COMPILER)' 'ifort'
  FORTRAN_COMPILER = ifort
else ifeq '$(MPIF90_COMPILER)' 'Cray'
  FORTRAN_COMPILER = crayftn
else
  $(error Unrecognised mpif90 compiler option: "$MPIF90_COMPILER")
endif

include $(LFRIC_BUILD)/fortran/$(FORTRAN_COMPILER).mk

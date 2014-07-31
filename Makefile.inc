##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Variables defined here will will be included in all Makefiles for this
# project.
# Expects ROOT to point to the root of the project tree.

ROOT ?= .
BUILD_DIR = $(ROOT)/build
TOOL_DIR = $(ROOT)/tools

ifdef PE_ENV
    # This may not be valid but we are going to assume that this variable
    # indicates a Cray XC-30.
    ifeq '$(PE_ENV)' 'CRAY'
        FC = crayftn
    else ifeq '$(PE_ENV)' 'INTEL'
        FC = ifort
    else ifeq '$(PE_ENV)' 'GNU'
        FC = gfortran
    else
        $(error Unrecognised Cray programming environment)
    endif
endif

ifndef FC
  OS = $(shell uname -s)
  MACHINE = $(shell uname -m)

  ifeq '$(OS)' 'Linux'
    ifeq '$(MACHINE)' 'x86_64'
      FC = ifort
    else
      $(error Unrecognised processor: $(MACHINE))
    endif
  else ifeq '$(OS)' 'AIX'
    FC = xlf
  else ifeq '$(OS)' 'Darwin'
    FC = gfortran
  else
    $(error Unrecognised operating system: $(OS))
  endif
endif

ifeq '$(FC)' 'ifort'
  $(info ** Chosen Intel Fortran compiler)

  FFLAGS_COMPILER     = -xhost
  FFLAGS_OPTIMISATION = -O0
  FFLAGS_DEBUG        = -g -traceback
  FFLAGS_WARNINGS     = -warn all -warn errors
  FFLAGS_RUNTIME      = -check all
  F_MOD_DESTINATION_ARG = -module $(OBJ_DIR)
  IFORT_VERSION = $(shell ifort -v 2>&1 | awk '{ print $$NF }' \
                   | awk -F . '{ printf "%03i%02i%02i", $$1, $$2, $$3 }')
  $(info ** Version $(IFORT_VERSION))

else ifeq '$(FC)' 'xlf'
  $(info ** Chosen IBM XL Fortran compiler)

  F_MOD_DESTINATION_ARG = -qmoddir=$(OBJ_DIR)

else ifeq '$(FC)' 'gfortran'
  $(info ** Chosen GNU Fortran compiler)

  FFLAGS_OPTIMISATION = -O0
  FFLAGS_DEBUG        = -g
  FFLAGS_WARNINGS     = -Wall -Werror
  FFLAGS_RUNTIME      = -fcheck=all
  F_MOD_DESTINATION_ARG = -J $(OBJ_DIR)
  GFORTRAN_VERSION = $(shell gfortran --version 2>&1 \
                       | sed -n 's/.*[[:space:]]\{1,\}\([[:digit:]]\{1,\}\.[[:digit:]]\{1,\}\.[[:digit:]]\{1,\}\)[[:space:]]\{,\}.*/\1/p' \
                       | awk -F . '{ printf "%02i%02i%02i", $$1, $$2, $$3 }')
  $(info ** Version $(GFORTRAN_VERSION))

else ifeq '$(FC)' 'nagfor'
  $(info ** Chosen NAG Fortran compiler)

  FFLAGS_OPTIMISATION = -O0
  FFLAGS_DEBUG        = -g -gline -colour
  FFLAGS_RUNTIME      =  -C=all -mtrace=all -nan
  F_MOD_DESTINATION_ARG = -mdir $(OBJ_DIR)

else ifeq '$(FC)' 'pgfortran'
  $(info ** Chosen Portland Fortran compiler)

  FFLAGS_OPTIMISATION = -O0
  FFLAGS_DEBUG        = -g -traceback
  FFLAGS_RUNTIME      = -Mbounds -Mchkptr -Mchkstk
  F_MOD_DESTINATION_ARG = -module $(OBJ_DIR)

else ifeq '$(FC)' 'crayftn'
  $(info ** Chosen Cray Fortran compiler)

  FFLAGS_COMPILER = -em
  FFLAGS_OPTIMISATION = -O0
  FFLAGS_DEBUG = -g
  FFLAGS_WARNINGS = -m 0
  F_MOD_DESTINATION_ARG = -J $(OBJ_DIR)

else
  $(error Unrecognised Fortran compiler $(FC))

endif

# Once we've set everything up we need to repair FC on Cray systems
ifdef PE_ENV
  FC = ftn
endif

FFLAGS += $(FFLAGS_COMPILER) $(FFLAGS_OPTIMISATION) $(FFLAGS_DEBUG) \
          $(FFLAGS_WARNINGS) $(FFLAGS_RUNTIME)

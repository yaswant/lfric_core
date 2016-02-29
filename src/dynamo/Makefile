##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

export IGNORE_DEPENDENCIES = netcdf MPI ESMF
export EXTERNAL_DYNAMIC_LIBRARIES = esmf netcdff netcdf
export EXTERNAL_STATIC_LIBRARIES  =

export DYNAMO_DEPENDENCY_DB ?= $(OBJ_DIR)/dependencies.db

export ROOT = ../..
export PSYCLONE_DIR = $(ROOT)/src/psyclone

include $(ROOT)/make/include.mk
include $(MAKE_DIR)/compiler.mk
include $(MAKE_DIR)/fortran/$(FORTRAN_COMPILER).mk
include $(MAKE_DIR)/mpi_link.mk
include make/fortran/$(FORTRAN_COMPILER).mk

OBJ_DIR = $(BUILD_DIR)/dynamo

LDFLAGS += $(FFLAGS_COMPILER)

ifdef CRAY_ENVIRONMENT
    $(info Building in Cray environment)
    # Cray's do not like linking dynamically so we make everything static.
    #
    override LINK := STATIC
endif

FFLAGS += $(FFLAGS_COMPILER) $(FFLAGS_WARNINGS)

ifeq ($(LINK),STATIC)
    export EXTERNAL_STATIC_LIBRARIES  := $(EXTERNAL_STATIC_LIBRARIES) \
                                         $(EXTERNAL_DYNAMIC_LIBRARIES)
    export EXTERNAL_DYNAMIC_LIBRARIES  = 
endif

ifeq ($(OPTIMISATION),NONE)
  FFLAGS += $(FFLAGS_NO_OPTIMISATION)
else ifeq ($(OPTIMISATION),SAFE)
  FFLAGS += $(FFLAGS_SAFE_OPTIMISATION)
else ifeq ($(OPTIMISATION),RISKY)
  FFLAGS += $(FFLAGS_RISKY_OPTIMISATION)
endif

ifeq ($(SYMBOLS),YES)
  FFLAGS += $(FFLAGS_DEBUG)
endif

ifeq ($(CHECKS),YES)
  FFLAGS += $(FFLAGS_INIT) $(FFLAGS_RUNTIME)
endif

export OBJ_DIR TOOL_DIR F_MOD_DESTINATION_ARG FFLAGS MPI_LD LDFLAGS DEPRULE_FLAGS

.PHONY: all
all: modules

.PHONY: modules
modules: applications
	$(MAKE) -f make/dynamo.mk modules

.PHONY: applications
applications:
	$(MAKE) -f make/dynamo.mk applications

.PHONY: unused-check
unused-check: | $(OBJ_DIR)
	$(MAKE) -f make/examine.mk unused-check

.PHONY: clean
clean:
	$(MAKE) -f make/dynamo.mk clean

$(OBJ_DIR):
	$(Q)mkdir -p $@
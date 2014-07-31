##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

ROOT = ../..
include $(ROOT)/Makefile.inc
FFLAGS := $(filter-out -Werror,$(FFLAGS)) # pFUnit isn't that well put together

TARGET_OBJ_DIR = $(BUILD_DIR)/dynamo

# This is not pretty, or good, but it is expedient. Until such times as I have
# written a better build system put special object files needed to link here.
LINK_SPECIAL = 

TEST_SUBDIRS = $(shell find * -maxdepth 1 -type d)

TEST_OBJ_DIR     = $(BUILD_DIR)/tests
TEST_OBJ_SUBDIRS = $(patsubst %,$(TEST_OBJ_DIR)/%,$(TEST_SUBDIRS))
# OBJ_DIR is used by Makefile.inc
OBJ_DIR          = $(TEST_OBJ_DIR)

TEST_SRC  = $(shell find . -name '*.pf')
TEST_OBJS = $(patsubst %.pf,$(TEST_OBJ_DIR)/%.o,$(TEST_SRC))

SRC  = $(wildcard *.[Ff]90)
OBJS = $(patsubst %.[Ff]90,$(TEST_OBJ_DIR)/%.o,$(SRC))

DRIVER_OBJ = $(TEST_OBJ_DIR)/driver.o
TARGET_OBJ = $(patsubst %_test.pf,$(TARGET_OBJ_DIR)/%.o,$(TEST_SRC))

DEP_FILE   = dependencies.mk

.PHONY: test
test: pfunit $(TEST_OBJ_DIR)/test
	$(TEST_OBJ_DIR)/test

include pfunit.mk

$(TEST_OBJ_DIR)/test: $(TEST_OBJS) $(OBJS) $(DRIVER_OBJ) $(TARGET_OBJ) \
                      $(LINK_SPECIAL)
	$(FC) -L$(PFUNIT_INSTALL)/lib -L/data/cr1/mhambley/modules/packages/gfortran/4.9.0/netcdf/4.3.2/lib \
        -o $@ $^ -lpfunit -lnetcdf -lnetcdff

$(DRIVER_OBJ): $(PFUNIT_INSTALL)/include/driver.F90 testSuites.inc
	$(FC) $(FFLAGS) -c -I$(PFUNIT_INSTALL)/mod -I. -o $@ $<

$(TEST_OBJ_DIR)/%.o: $(TEST_OBJ_DIR)/%.F90
	@echo Compiling $<
	$(FC) $(FFLAGS) -c $(F_MOD_DESTINATION_ARG) \
	    -I$(OBJ_DIR) -I$(PFUNIT_INSTALL)/mod \
	      -I$(TARGET_OBJ_DIR) -o $@ $<

.PRECIOUS: $(TEST_OBJ_DIR)/%.F90
$(TEST_OBJ_DIR)/%.F90: %.pf $(TEST_OBJ_DIR) $(TEST_OBJ_SUBDIRS)
	@echo Generating $@
	$(PFUNIT_INSTALL)/bin/pFUnitParser.py $< $@

$(TEST_OBJ_DIR) $(TEST_OBJ_SUBDIRS):
	@echo "Creating $@"
	@mkdir -p $@

.PHONY: clean
clean: cleanpfunit
	-rm -rf $(OBJ_DIR)

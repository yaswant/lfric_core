##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

ALGORITHM_PATH = algorithm
PSY_PATH       = psy
KERNEL_PATH    = kernel

AUTO_ALGORITHM_PATH = $(OBJ_DIR)/algorithm
AUTO_PSY_PATH       = $(OBJ_DIR)/psy

PSYCLONE = python2.7 $(PSYCLONE_DIR)/src/generator.py

include $(ROOT)/make/include.mk

AUTO_FILES   := $(patsubst $(ALGORITHM_PATH)/%.x90, \
                           $(AUTO_PSY_PATH)/psy_%.f90, \
                           $(wildcard $(ALGORITHM_PATH)/*.x90) )
MANUAL_FILES := $(wildcard $(PSY_PATH)/*.[Ff]90 )
AUTO_FILES   := $(filter-out $(patsubst $(PSY_PATH)/%, \
                                        $(AUTO_PSY_PATH)/%, \
                                        $(MANUAL_FILES) ), \
                             $(AUTO_FILES) )
MANUAL_FILES := $(patsubst $(PSY_PATH)/psy_%, \
                           $(AUTO_ALGORITHM_PATH)/%, \
                           $(MANUAL_FILES) )

# The appending operator (+=) is not used below as Cylc adds a space on the
# end of the variable, thereby ruining anyone elses chance to append.
# Instead we prepend manually.
#
export PYTHONPATH := $(PSYCLONE_DIR)/f2py_93:$(PSYCLONE_DIR)/src:$(PYTHONPATH)

.PHONY: all
all: $(AUTO_FILES) $(MANUAL_FILES)

$(AUTO_ALGORITHM_PATH)/%.f90:$(ALGORITHM_PATH)/%.x90 | $(AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)PSycloning$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -d $(KERNEL_PATH) \
	            -oalg $(patsubst $(ALGORITHM_PATH)/%.x90, $(AUTO_ALGORITHM_PATH)/%.f90, $< ) \
	            $<

$(AUTO_PSY_PATH)/psy_%.f90: $(ALGORITHM_PATH)/%.x90 | $(AUTO_PSY_PATH) $(AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)PSycloning$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -d $(KERNEL_PATH) \
	            -opsy $(patsubst $(ALGORITHM_PATH)/%.x90, \
	                             $(AUTO_PSY_PATH)/psy_%.f90, $< ) \
	            -oalg $(patsubst $(ALGORITHM_PATH)/%.x90, \
	                             $(AUTO_ALGORITHM_PATH)/%.f90, $< ) \
	            $<

$(AUTO_PSY_PATH) $(AUTO_ALGORITHM_PATH):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

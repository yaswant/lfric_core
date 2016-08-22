##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

PSYCLONE ?= python2.7 $(PSYCLONE_DIR)/src/generator.py

PSY_ALGORITHM_PATH = algorithm
PSY_PSY_PATH       = psy
PSY_KERNEL_PATH    = kernel

PSY_OPTIMISATION_PATH = $(ROOT)/optimisation

PSY_AUTO_ALGORITHM_PATH = $(OBJ_DIR)/algorithm
PSY_AUTO_PSY_PATH       = $(OBJ_DIR)/psy

PSY_AUTO_FILES   := $(patsubst $(PSY_ALGORITHM_PATH)/%.x90, \
                               $(PSY_AUTO_PSY_PATH)/psy_%.f90, \
                               $(wildcard $(PSY_ALGORITHM_PATH)/*.x90) )
PSY_MANUAL_FILES := $(wildcard $(PSY_PSY_PATH)/*.[Ff]90 )
PSY_AUTO_FILES   := $(filter-out $(patsubst $(PSY_PSY_PATH)/%, \
                                        $(PSY_AUTO_PSY_PATH)/%, \
                                        $(PSY_MANUAL_FILES) ), \
                             $(PSY_AUTO_FILES) )
PSY_MANUAL_FILES := $(patsubst $(PSY_PSY_PATH)/psy_%, \
                           $(PSY_AUTO_ALGORITHM_PATH)/%, \
                           $(PSY_MANUAL_FILES) )

# The appending operator (+=) is not used below as Cylc adds a space on the
# end of the variable, thereby ruining anyone elses chance to append.
# Instead we prepend manually.
#
export PYTHONPATH := $(PSYCLONE_DIR)/f2py_93:$(PSYCLONE_DIR)/src:$(PYTHONPATH)

.PHONY: generate-psykal
generate-psykal: $(PSY_AUTO_FILES) $(PSY_MANUAL_FILES) $(PSY_ALGORITHM_FILES)

$(PSY_AUTO_ALGORITHM_PATH)/%.f90: $(PSY_ALGORITHM_PATH)/%.x90 \
                                  | $(PSY_AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)Overridden PSyclone$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -l -d $(PSY_KERNEL_PATH) \
	            -oalg $@ \
	            $<

$(PSY_AUTO_PSY_PATH)/psy_%.f90: $(PSY_ALGORITHM_PATH)/%.x90  \
                        $(PSY_OPTIMISATION_PATH)/$(DYNAMO_OPTIMISATION_PROFILE)/%.py \
                                | $(PSY_AUTO_PSY_PATH)       \
                                  $(PSY_AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)Full PSyclone, local optimisations$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -l -d $(PSY_KERNEL_PATH)                     \
	            -s $(PSY_OPTIMISATION_PATH)/$(DYNAMO_OPTIMISATION_PROFILE)/$*.py \
	            -opsy $@                                                 \
	            -oalg $(patsubst $(PSY_ALGORITHM_PATH)/%.x90, \
	                             $(PSY_AUTO_ALGORITHM_PATH)/%.f90, $< )  \
	            $<

$(PSY_AUTO_PSY_PATH)/psy_%.f90: $(PSY_ALGORITHM_PATH)/%.x90  \
                    $(PSY_OPTIMISATION_PATH)/$(DYNAMO_OPTIMISATION_PROFILE)/global.py \
                                | $(PSY_AUTO_PSY_PATH)       \
                                  $(PSY_AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)Full PSyclone, global optimisations$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -l -d $(PSY_KERNEL_PATH)                      \
	         -s $(PSY_OPTIMISATION_PATH)/$(DYNAMO_OPTIMISATION_PROFILE)/global.py \
	            -opsy $@                                                  \
	            -oalg $(patsubst $(PSY_ALGORITHM_PATH)/%.x90, \
	                             $(PSY_AUTO_ALGORITHM_PATH)/%.f90, $< )   \
	            $<

$(PSY_AUTO_PSY_PATH)/psy_%.f90: $(PSY_ALGORITHM_PATH)/%.x90 \
                                | $(PSY_AUTO_PSY_PATH)      \
                                  $(PSY_AUTO_ALGORITHM_PATH)
	@echo -e $(VT_BOLD)Full PSyclone$(VT_RESET) $<
	$(PSYCLONE) -api dynamo0.3 -l -d $(PSY_KERNEL_PATH) \
	            -opsy $@ \
	            -oalg $(patsubst $(PSY_ALGORITHM_PATH)/%.x90, \
	                             $(PSY_AUTO_ALGORITHM_PATH)/%.f90, $< ) \
	            $<

$(PSY_AUTO_PSY_PATH) $(PSY_AUTO_ALGORITHM_PATH):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

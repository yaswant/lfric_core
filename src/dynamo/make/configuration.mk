##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

CONFIGURATOR ?= $(TOOL_DIR)/GenerateNamelist
CONFIGULOADER ?= $(TOOL_DIR)/GenerateLoader

CONFIGURATION_PATH = configuration
AUTO_CONFIGURATION_PATH = $(OBJ_DIR)/configuration

CONFIGURATION_SOURCE_FILES := $(patsubst $(CONFIGURATION_PATH)/%.nld,\
                           $(AUTO_CONFIGURATION_PATH)/%_config_mod.f90,\
                           $(wildcard $(CONFIGURATION_PATH)/*.nld))

CONFIGURATION_NAMELISTS := $(patsubst $(CONFIGURATION_PATH)/%.nld, %, \
                                      $(wildcard $(CONFIGURATION_PATH)/*.nld))

CONFIGURATION_READER_SOURCE = $(AUTO_CONFIGURATION_PATH)/configuration_mod.f90

.PHONY: generate-configuration
generate-configuration: $(CONFIGURATION_SOURCE_FILES) \
                        $(CONFIGURATION_READER_SOURCE)

$(AUTO_CONFIGURATION_PATH)/%_config_mod.f90: $(CONFIGURATION_PATH)/%.nld\
                                             | $(AUTO_CONFIGURATION_PATH)
	@echo -e $(VT_BOLD)Configurator$(VT_RESET) $<
	$(CONFIGURATOR) -directory $(AUTO_CONFIGURATION_PATH) $<

$(CONFIGURATION_READER_SOURCE): $(CONFIGURATION_SOURCE_FILES)
	@echo -e $(VT_BOLD)Configuration loader$(VT_RESET) $@
	$(CONFIGULOADER) $@ $(CONFIGURATION_NAMELISTS)

$(AUTO_CONFIGURATION_PATH):
	@echo Creating $@
	@mkdir -p $@
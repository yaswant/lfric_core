##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

DESIGN_DOCUMENTS = formulation/dynamo_formulation.lyx glossary/Glossary.lyx
DESIGN_OUTPUTS = $(patsubst %.lyx,%.pdf,$(DESIGN_DOCUMENTS))

DOXYGEN ?= doxygen

TOP := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
$(info ** TEXINPUTS = $(TOP))

.PHONY: most
most: api design

.PHONY: api
api:
	echo Building API documentation
	$(DOXYGEN) Doxyfile

.PHONY: design
design: $(DESIGN_OUTPUTS)

%.pdf: %.lyx
	echo Building $@
	env TEXINPUTS=$(TOP)/tex//: lyx -e pdf $<

.PHONY: gource
gource:
	$(MAKE) -C visualisation

.PHONY: clean
clean:
	-rm -rf docs/api/*
	-rm $(DESIGN_OUTPUTS)


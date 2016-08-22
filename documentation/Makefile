##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

DOXYGEN  ?= doxygen

DOCUMENT_SOURCES = $(wildcard */*.latex)
DOCUMENT_DIRECTORIES = $(patsubst %/,%,$(dir $(DOCUMENT_SOURCES)))
API_COLLECTIONS = $(patsubst Doxyfile.%,api/%,$(wildcard Doxyfile.*))

DYNAMO_COMMON_FIGURES := $(shell pwd)/common-figures
export DYNAMO_COMMON_FIGURES

.PHONY: most
most: design-documents api-documents

.PHONY: design-documents
design-documents: $(DOCUMENT_DIRECTORIES) uml

design-documents.tar: design-documents
	gtar --mode=644 -cf $@ uml/*.pdf
	for filename in $(patsubst %.latex,%.pdf,$(DOCUMENT_SOURCES)); do \
	    gtar --mode=644 -C `dirname $$filename` -rf $@  `basename $$filename`; \
	done

.PHONY: api-documents
api-documents: $(API_COLLECTIONS)

.PHONY: api/%
api/%: Doxyfile.% | api
	echo Building $* API documentation
	$(DOXYGEN) $<

api:
	echo Creating directory "$@"
	mkdir $@

.PHONY: $(DOCUMENT_DIRECTORIES)
$(DOCUMENT_DIRECTORIES):
	$(MAKE) -f ../make/latex.mk -C $@

.PHONY: uml
uml:
	$(MAKE) -C uml

.PHONY: gource
gource:
	$(MAKE) -C visualisation

clean: $(addprefix clean_,$(DOCUMENT_DIRECTORIES)) clean_uml
	-rm -rf api/*

cleam_uml: ALWAYS
	$(MAKE) -C uml clean

clean_%: ALWAYS
	$(MAKE) -f ../make/latex.mk -C $* clean

.PHONY: ALWAYS
ALWAYS:


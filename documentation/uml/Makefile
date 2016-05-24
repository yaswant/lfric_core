DIAGRAMS = $(patsubst %.puml,%.pdf,$(wildcard *.puml))

.PHONY: all
all: $(DIAGRAMS)

.SECONDEXPANSION:

%.pdf: %.svg
	inkscape $< --export-pdf=$@

%.svg: %.puml $$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$*.puml)
	plantuml -tsvg $<

.PHONY: clean
clean:
	-rm -f *.svg *.pdf

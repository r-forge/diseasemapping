ifndef R_HOME
	R_HOME:=/usr
endif

PDFS:= mapmisc.pdf mapLayers.pdf
HTML:= mapLayers.html

Rargs? = ""


all: $(PDFS) 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref

%.tex: %.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" $(Rargs)

%.pdf: %.tex
	xelatex $<;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')" 

%.html: %.tex
	htlatex $<;

web:	mapLayers.html
	cp -r figure ../../../www/openmap;
	cp mapLayers.html ../../../www/openmap/index.html;
	cp mapLayers.css ../../../www/openmap/


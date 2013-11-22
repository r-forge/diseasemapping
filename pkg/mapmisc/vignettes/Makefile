ifndef R_HOME
	R_HOME:=/usr
endif

PDFS:= mapmisc.pdf mapLayers.pdf
HTML:= mapLayers.html

all: $(PDFS) 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/

%.tex: %.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" 

%.pdf: %.tex
	xelatex $<;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')" 

%.html: %.rnw
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2html('$<.rnw')"

web:	mapLayers.html
	cp -r figure/someTiles*.png ../../www/openmap/figure/;
	cp mapLayers.html ../../www/openmap/index.html;
	cp mapLayers.css ../../www/openmap/


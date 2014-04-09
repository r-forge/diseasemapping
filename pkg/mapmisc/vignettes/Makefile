ifndef R_HOME
	R_HOME:=/usr
endif

PDFS:= mapmisc.pdf mapLayers.pdf

Rargs? = ""


all: $(PDFS) 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref

%.tex: %.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" $(Rargs)

%.pdf: %.tex
	xelatex $<;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')" 


web:	mapLayers.tex
	cp mapLayers.tex mapWeb.tex
	htlatex mapWeb.tex
	cp -r figure ../../../www/openmap;
	cp mapWeb.html ../../../www/openmap/index.html;
	cp mapWeb.css ../../../www/openmap/


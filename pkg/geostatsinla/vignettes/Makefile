ifndef R_HOME
	R_HOME:=/usr
endif

PDFS:= lgcp.pdf

all: $(PDFS) 

rargs? = ""


clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.lg *.RData

%.tex: %.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" $(rargs)

%.pdf: %.tex
	xelatex $<;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')" 


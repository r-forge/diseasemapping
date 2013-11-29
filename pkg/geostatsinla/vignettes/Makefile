ifndef R_HOME
	R_HOME:=/usr
endif

PDFS:= lgcp.pdf

all: $(PDFS) 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi

%.tex: %.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" 

%.pdf: %.tex
	xelatex $<;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('$@', gs_quality='ebook')" 


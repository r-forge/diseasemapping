ifndef R_HOME
	R_HOME:=/usr
endif


all: lgcp.pdf

Rargs? = "reCache smallMemory"

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.lg *.RData

lgcp.tex: lgcp.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('lgcp.rnw')" $(Rargs)

varioustex: lgcp.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" reRun reCache smallMemory;
	xelatex lgcp.tex;
	mv lgcp.pdf lgcpReRunReCachesmallMemory.pdf;
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" reRun smallMemory;
	xelatex lgcp.tex;
	mv lgcp.pdf lgcpReRunsmallMemory.pdf;
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" smallMemory;
	xelatex lgcp.tex;
	mv lgcp.pdf lgcpsmallMemory.pdf;
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" reCache smallMemory;
	xelatex lgcp.tex;
	mv lgcp.pdf lgcpreCachesmallMemory.pdf;
	rm lgcp.tex
 
resave: lgcp.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('$<')" reRun reCache reSave;

lgcp.pdf: lgcp.tex
	xelatex lgcp;
	$(R_HOME)/bin/Rscript -e "tools::compactPDF('lgcp.pdf', gs_quality='ebook')" 


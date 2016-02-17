all: bym.pdf

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz

bym.Rnw: bym.R
	'$(R_HOME)/bin/Rscript' -e "knitr::spin('bym.R', knit=FALSE, format='Rnw')" $(Rargs)


bym.tex: bym.Rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('bym.Rnw')" $(Rargs)

bym.pdf: bym.tex
	sed '2 a \\\\usepackage{subcaption}' bym.tex > bym2.tex;
	xelatex bym2;
	xelatex bym2;
	cp bym2.pdf bym.pdf;
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('bym.pdf', gs_quality='ebook')" 



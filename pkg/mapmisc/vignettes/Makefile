all: mapLayers.html north.html mapmisc.pdf 

clean:
	rm -rf *.md *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.pdf *.html

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')" 

north.md: north.Rmd
	R --slave -e "knitr::knit('north.Rmd')" fact=2
	
%.html: %.md
	pandoc --smart --self-contained --to=html5 --output=$@ $<	
all: mapLayers.html north.html mapmisc.pdf 

clean:
	rm -rf *.md *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.pdf *.html

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')" 

%.html: %.md
	pandoc --standalone --self-contained --smart --output=$@ $<

%.tex: %.Rnw
	R --slave -e "knitr::knit('$<')" $(Rargs)

%.md: %.Rmd
	R --slave -e "knitr::knit('$<')" $(Rargs)

web: mapmisc.pdf north.html mapLayers.html
	cp mapmisc.pdf ../../../www/mapmisc.pdf;
	cp north.html ../../../www/north.html;
	cp mapLayers.html ../../../www/openmap.html

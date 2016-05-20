all: mapLayers.html mapmisc.pdf 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz


mapLayers.html: mapLayers.Rmd
	'$(R_HOME)/bin/Rscript' -e "knitr::knit2html('mapLayers.Rmd')" $(Rargs)

mapmisc.tex: mapmisc.rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('mapmisc.rnw')" $(Rargs)

north.html: north.Rmd
	'$(R_HOME)/bin/Rscript' -e "knitr::knit2html('north.Rmd')" $(Rargs)

mapmisc.pdf: mapmisc.tex
	xelatex mapmisc;
	xelatex mapmisc
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('mapmisc.pdf', gs_quality='ebook')" 

web:
	cp mapmisc.pdf ../../../www/mapmisc.pdf;
	cp north.html ../../../www/north.html;
	cp mapLayers.html ../../../www/openmap.html

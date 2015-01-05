all: mapLayers.pdf mapmisc.pdf 

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz


mapLayers.tex: mapLayers.rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('mapLayers.rnw')" $(Rargs)

mapmisc.tex: mapmisc.rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('mapmisc.rnw')" $(Rargs)

mapmisc.pdf: mapmisc.tex
	xelatex mapmisc;
	xelatex mapmisc
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('mapmisc.pdf', gs_quality='ebook')" 

mapLayers.pdf: mapLayers.tex
	xelatex mapLayers;
	xelatex mapLayers
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('mapLayers.pdf', gs_quality='ebook')" 

web:	mapLayers.tex
	cp mapLayers.tex mapWeb.tex
	htlatex mapWeb.tex
	cp -r figure ../../../www/openmap;
	cp mapWeb.html ../../../www/openmap/index.html;
	cp mapWeb.css ../../../www/openmap/

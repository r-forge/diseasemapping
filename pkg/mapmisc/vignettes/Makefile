all: mapLayers.pdf mapLayers.html mapmisc.pdf

mapLayers.pdf: mapLayers.tex
	xelatex mapLayers.tex;
	xelatex mapLayers.tex

mapLayers.html: mapLayers.tex
	htlatex mapLayers.tex;
	htlatex mapLayers.tex

mapLayers.tex: mapLayers.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('mapLayers.rnw')"

mapmisc.tex: mapmisc.rnw
	$(R_HOME)/bin/Rscript -e "knitr::knit('mapmisc.rnw')"

mapmisc.pdf: mapmisc.tex
	xelatex mapmisc.tex;
	xelatex mapmisc.tex

web:	mapLayers.html
	cp -r figure ../../../www/openmap/;
	cp mapLayers.html ../../../www/openmap/index.html;
	cp mapLayers.css ../../../www/openmap/

clean:	mapmisc.pdf mapLayers.html mapLayers.pdf	
	rm *.dvi *.toc *.log *.aux *.tex *.lg *.idv *.4* *.tmp *.xref

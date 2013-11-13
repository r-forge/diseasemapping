PDFS= mapmisc.pdf mapLayers.pdf
HTML = mapLayers.html

all: $(PDFS) 

clean:
        rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/

%.pdf: %.Rnw
        $(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2pdf('$*.Rnw') else tools::texi2pdf('$*.tex')"

%.html: %.Rmd
        $(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2html('$*.Rmd')"

web:	mapLayers.html
	cp -r figure/someTiles*.png ../../../www/openmap/figure/;
	cp mapLayers.html ../../../www/openmap/index.html;
	cp mapLayers.css ../../../www/openmap/


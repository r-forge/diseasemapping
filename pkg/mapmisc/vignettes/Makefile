all: mapLayers.pdf mapLayers.html

mapLayers.pdf: mapLayers.tex
	pdflatex mapLayers.tex

mapLayers.html: mapLayers.tex
	htlatex mapLayers.tex

mapmisc.pdf: mapmisc.tex
	pdflatex mapmisc.tex

web:	mapLayers.html
	cp -r figure ../../../www/openmap/;
	cp mapLayers.html ../../../www/openmap/index.html;
	cp mapLayers.css ../../../www/openmap/


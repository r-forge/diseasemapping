VPATH = ../pkg/mapmisc/inst/extdoc:../pkg/mapmisc/vignettes:../pkg/mapmisc/tests/


finalDocs = world.html mapLayers.html north.html mapmisc.pdf openmap.html worldMap.html
testDocs = test_perspective.html test_north.html test_openmap.html test_worldMap.html

all: $(finalDocs)

test: $(testDocs)

openmap.Rmd: openmap.R
	R --slave -e "knitr::spin('$^', knit=FALSE)";
	mv ../pkg/mapmisc/tests/openmap.Rmd .

worldMap.Rmd: worldMap.R
	R --slave -e "knitr::spin('$^', knit=FALSE)";
	mv ../pkg/mapmisc/tests/worldMap.Rmd .

clean:
	rm -rf *.md *.tex *.bcf *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.pdf *.html

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')" 

%.tex: %.md
	pandoc --standalone --smart --biblatex --from=markdown --output=$@ $<

%.md: %.Rmd
	R --slave -e "testing=FALSE;fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=FALSE;knitr::knit('$^', '$@')"

test_%.md: %.Rmd
	R --slave -e "testing=TRUE;fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=TRUE;knitr::knit('$^', '$@')"

%.tex: %.Rnw
	R --slave -e "fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=FALSE;knitr::knit('$<', '$@')" 
	
%.html: %.md
	pandoc --smart --standalone --to=html5 --css=http://diseasemapping.r-forge.r-project.org/webpage.css --output=$@ $<	


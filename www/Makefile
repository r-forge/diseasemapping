VPATH = ../pkg/geostatsp/vignettes

finalDocs = glgm.pdf north.html mapmisc.pdf openmap.html worldMap.html mapLayers.html perspective.html
testDocs = test_perspective.html test_north.html test_openmap.html test_worldMap.html test_mapmisc.pdf


all: $(finalDocs)
test: $(testDocs)

testExistingDocs: test_existing_cache_perspective.html test_existing_cache_north.html test_existing_cache_openmap.html test_existing_cache_worldMap.html test_existing_cache_mapmisc.pdf

%.Rmd: ../pkg/mapmisc/inst/extdoc/%.Rmd
	cp $< $@

%.Rmd: ../pkg/mapmisc/vignettes/%.Rmd
	cp $< $@

%.Rnw: ../pkg/mapmisc/vignettes/%.Rnw
	cp $< $@


%.R: ../pkg/mapmisc/tests/%.R
	mv $< .;

%.Rmd: %.R
	R --slave -e "knitr::spin('$<', knit=FALSE)";

clean:
	rm -rf *.md *.tex *.bcf *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm -r *.pdf *.html mapmiscCache

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')" 

%.tex: %.md
	pandoc --standalone --smart --biblatex --from=markdown --output=$@ $<

%.md: %.Rmd
	R --slave -e "testing=FALSE;fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=FALSE;knitr::knit('$<', '$@')"

test_%.md: %.Rmd %.html
	R --slave -e "testing=TRUE;fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=TRUE;knitr::knit('$<', '$@')"

test_%.tex: %.Rnw %.pdf
	R --slave -e "testing=TRUE;fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=TRUE;knitr::knit('$<', '$@')"

test_existing_cache_%.md: %.Rmd
	R --slave -e "testing=TRUE;fact=2;mapmiscCacheReadOnly=TRUE;knitr::knit('$<', '$@')"

test_existing_cache_%.tex: %.Rnw
	R --slave -e "testing=TRUE;fact=2;mapmiscCacheReadOnly=TRUE;knitr::knit('$<', '$@')"

%.tex: %.Rnw
	R --slave -e "fact=2;mapmiscCachePath=file.path(getwd(), 'mapmiscCache');mapmiscCacheReadOnly=FALSE;knitr::knit('$<', '$@')" 
	
%.html: %.md
	pandoc --smart --standalone --to=html5 --css=http://diseasemapping.r-forge.r-project.org/webpage.css --output=$@ $<	



# to do: mapLayers.Rmd is lost
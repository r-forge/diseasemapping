build: build-compact

build-compact:
	cd .. && R CMD build --compact-vignettes=both mapmiscTerra

compact-vignettes:
	Rscript -e "tools::compactPDF(paths='inst/doc', gs_quality='ebook')"

check: build-compact
	R CMD check --as-cran ../*.tar.gz

clean:
	rm -f ../*.tar.gz

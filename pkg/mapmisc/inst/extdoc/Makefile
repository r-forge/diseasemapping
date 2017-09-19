all: world.pdf 

clean:
	rm -rf *.md *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.pdf *.html

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<

%.md: %.Rmd
	R --slave -e "knitr::knit('$<')"
	
%.html: %.md
	pandoc --smart --self-contained --filter=pandoc-crossref --to=html5 --output=$@ $<

%.tex: %.md
	pandoc --standalone --smart --filter=pandoc-crossref --from=markdown --to=latex $< | sed s/\\\\usepackage{subfig}// > $@
	
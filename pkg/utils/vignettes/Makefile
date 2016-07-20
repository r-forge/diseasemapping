all: slides.pdf

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz

%.md: %.Rmd
	R -e "knitr::knit('$<', encoding='UTF-8')" 

%.tex: %.md
	pandoc --standalone --smart --to=beamer --output=$@ $<

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<


